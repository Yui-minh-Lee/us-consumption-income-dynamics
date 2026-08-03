%% Step 5: 四变量-双周期-固定参数SSM
% 结构: 
%   y_t = tau_y + psi_t
%   c_t = tau_y + lambda_t * psi_t + gamma * Tr_t
%   u_gap_t = theta * psi_t  (直接观测失业率缺口)
%   fiscal_cycle
% 目的: 利用失业率作为"锚"，强行纠正 2008 年的周期识别错误
% Plus:
%   1. 使用 "市场收入" (W875RX1) 识别纯粹的经济周期。
%   2. 将 "转移支付" (Trans_t) 独立出来解释消费，避免模型将财政刺激误判为趋势收入。
%   3. 现在y是市场收入 W875RX1。
%   4. 失业率趋势 (NAIRU) 在模型外通过 HP 滤波移除，不再作为内部状态估计。
%   5. 双周期模型。
%   6. 去掉时变敏感性

clear; clc; close all;

project_root = fileparts(fileparts(mfilename('fullpath')));

% =========================================================================
% 1. 数据准备 (复用Step2)
% =========================================================================
load(fullfile(project_root, 'data', 'processed', 'prepared_data.mat'));

% ---------- A. Hamilton滤波清洗：去趋势、整理协整关系 ----------
lead_len = 24;
lag_len = 12;
% 首部会有缺失, 截断起始NaN
data_time = data.Date(lead_len + lag_len : end);
T_len = length(data_time);

fprintf('正在使用 Hamilton 滤波进行去趋势处理 (LeadLength = %d, LagLength = %d)...\n', lead_len, lag_len);

% 1. 收入，不用去趋势 
[y_trend_h, y_cycle_h] = hfilter(ln_y, 'LeadLength', lead_len, 'LagLength', lag_len);
ln_y_adj = y_cycle_h(lead_len + lag_len : end); % 截断头部缺失值

% 2. 消费，去APC制度性漂移
apc_ratio = ln_c - ln_y;
[apc_trend_h, apc_cycle_h] = hfilter(apc_ratio, 'LeadLength', lead_len, 'LagLength', lag_len);
c_cycle_h = y_cycle_h + apc_cycle_h;
ln_c_adj = c_cycle_h(lead_len + lag_len : end);


% 3. 财政因子去趋势 (只保留刺激脉冲)
[trans_trend_h, trans_cycle_h] = hfilter(trans_factor, 'LeadLength', lead_len, 'LagLength', lag_len);
trans_factor_adj = trans_cycle_h(lead_len + lag_len : end);
fprintf('数据清洗完成。财政因子波动范围: %.4f ~ %.4f\n', min(trans_factor_adj), max(trans_factor_adj));


% 4. 失业率去趋势 HF滤波
unRate = data.Unemploy_rate / 100;
[u_trend_h, u_gap_h] = hfilter(unRate, 'LeadLength', lead_len, 'LagLength', lag_len);
unRate = unRate(lead_len + lag_len : end);
u_trend_h = u_trend_h(lead_len + lag_len : end);
u_gap_h = u_gap_h(lead_len + lag_len : end);


Data = [ln_y_adj, ln_c_adj, u_gap_h];


% ---------- B. 准备外部状态变量 S  --------
% S 必须是 "t-1" 时刻的信息，或者是当期外生的金融指标
S_raw = data.Unemploy_rate(lead_len + lag_len : end);
risk_threshold = mean(S_raw) + 1 * std(S_raw);  % threshold：1 std
S_normalized = double(S_raw > risk_threshold); % dummy array
S_normalized = [0; S_normalized(1:end-1)];       % lag 1

fprintf('状态变量 S_{t-1} 已构建并标准化 (Mean=%.2f, Std=%.2f)\n', mean(S_normalized), std(S_normalized));
fprintf('危机状态 S_{t-1} 占比: %.1f%%\n', mean(S_normalized)*100);


% =========================================================================
% 2. 继承 Step 2 的参数作为初始值
% =========================================================================

val_rho_long = 0.9042 ;   % long term damping
val_period_long = 80.0;  % long cycle period (month)
val_scyc_long = 0.00556;  % long cycle vol
val_lambda_long = 0.6063;  % long term lensitivity (const part)
val_okun_long = -0.7;     % 奥肯系数 (预期为负)

val_rho_short = 0.8 ;   % short term damping
val_period_short = 20;  % short cycle period (month)
val_scyc_short = 0.0050;  % short cycle vol
val_lambda_short = 0.5;  % short term sensitivity (const part)
val_okun_short = -0.3;     % 奥肯系数 (预期为负)

val_suy      = 0.00001;   % income noise
val_suc      = 0.0020;   % consump noise
%val_suu       = 0.0020;      % 失业率测量误差

% --- 逆变换：带边界保护---- 
% 将实际经济参数转回无约束空间，作为 fminunc 的起点

% 防止 log(0) 或 log(negative)
clip = @(x, min_v, max_v) max(min_v + 1e-6, min(max_v - 1e-6, x));
inv_sigmoid = @(x, min_v, max_v) -log((max_v - min_v)/(x - min_v) - 1);


% --- 长周期 ---
% 1. rho_long
rho_L_min = 0.5;
rho_L_max = 0.995;
p1 = inv_sigmoid(clip(val_rho_long, rho_L_min, rho_L_max), rho_L_min, rho_L_max);
% 2. freq_long (限制在 48-96 个月)
freq_L_min = 2*pi/96;
freq_L_max = 2*pi/60;
p2 = inv_sigmoid(clip(2*pi/val_period_long, freq_L_min, freq_L_max), freq_L_min, freq_L_max);
% 3. scyc_long
p3 = log(val_scyc_long);
% 4. lambda_long (预期接近 1)
p4 = val_lambda_long; 
% 5. theta_long (Okun, 预期为负)
p5 = log(-val_okun_long); % 存正数，Mapping取负

% --- 短周期---
% 6. rho_short (通常阻尼更小，消失得更快)
rho_S_min = 0.1;
rho_S_max = 0.95;
p6 = inv_sigmoid(clip(val_rho_short, rho_S_min, rho_S_max), rho_S_min, rho_S_max);
% 7. freq_short (限制在 12-32 个月)
freq_S_min = 2*pi/24;
freq_S_max = 2*pi/8;
p7 = inv_sigmoid(clip(2*pi/val_period_short, freq_S_min, freq_S_max), freq_S_min, freq_S_max);
% 8. scyc_short
p8 = log(val_scyc_short);
% 9. lambda_short
p9 = val_lambda_short; 
% 10. theta_short (短周期对失业率影响可能较小)
p10 = val_okun_short;

% 11. s_uc (exp)
p11 = log(val_suc);

% 12. U Noise (Log): > 0 constraint
%p12 = log(val_suu);


params0 = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11];


% =========================================================================
% 3. SSM 模型
% =========================================================================
Mdl_Double = ssm(@(p) my_double_cycle_mapping(p, trans_factor_adj, T_len));

% =========================================================================
% 4. 模型估计 (改进版：fminsearch + fminunc 混合策略)
% =========================================================================
fprintf('开始 Step 4 双时变参数估计 (TVP-SSM)...\n');
% --- 阶段 A: 使用 fminsearch 单纯形法进行搜索 ---
fprintf('\n--- Phase 1: Pre-optimization with fminsearch (Simplex) ---\n');
ObjFun = @(p) get_neg_logL(p, Data, trans_factor_adj, T_len);  % 定义目标函数句柄 (负对数似然)
options_search = optimset('MaxFunEvals',50000,'MaxIter',5000);  % 设置 fminsearch 选项
t_start = tic; % 记录时间
[params_search, fval_search] = fminsearch(ObjFun, params0, options_search); % 运行 fminsearch
fprintf('fminsearch 完成! 耗时: %.1f 秒, Min NLogL: %.4f\n', toc(t_start), fval_search);

% --- 阶段 B: 将 fminsearch 的结果作为 fminunc 的起点 ---
fprintf('\n--- Phase 2: Polishing with estimate (fminunc) ---\n');
fprintf('使用 fminsearch 找到的最优解作为新起点，开始计算 Hessian 和标准误...\n');

params0 = params_search; % 更新初始值
options = optimoptions('fminunc', 'MaxFunctionEvaluations', 50000, ...
    'MaxIterations', 5000, 'Display', 'iter', ...
    'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6);

try
    [EstMdl, estParams, EstCov, logL] = ...
        estimate(Mdl_Double, Data, params0, 'Options', options);
    fprintf('Step 4 最终估计完成！LogL: %.4f\n', logL);
catch ME
    fprintf('估计失败: %s\n', ME.message);
    return;
end

% =========================================================================
%% 5. 结果解析与绘图
% =========================================================================
% 提取状态
[SmoothedStates, LogL, Output] = smooth(EstMdl, Data);
SmoothedStatesCov = Output.SmoothedStatesCov;

% 状态顺序: [psi_L; psi*_L; psi_S; psi*_S; gamma]
psi_long_est  = SmoothedStates(:, 1); % 长周期
psi_short_est = SmoothedStates(:, 3); % 短周期
gamma_est_seq = SmoothedStates(:, 5); % 财政乘数
gamma_u_est_seq = SmoothedStates(:, 6); % 财政乘数

% 长周期参数
lambda_long = estParams(4);
theta_long  = -exp(estParams(5)); % 奥肯系数
rho_long    = rho_L_min + (rho_L_max - rho_L_min) / (1 + exp(-estParams(1)));
period_long = 2*pi / (freq_L_min + (freq_L_max - freq_L_min)/(1+exp(-estParams(2))));
scyc_long = exp(estParams(3));

% 短周期参数
lambda_short = estParams(9);
theta_short  = estParams(10); % 奥肯系数
rho_short    = rho_S_min + (rho_S_max - rho_S_min) / (1 + exp(-estParams(6)));
period_short = 2*pi / (freq_S_min + (freq_S_max - freq_S_min)/(1+exp(-estParams(7))));
scyc_short = exp(estParams(8));

% 财政乘数/敏感度
gamma_val  = mean(gamma_est_seq);
gamma_u_val  = mean(gamma_u_est_seq);

% 误差项方差
uc_val = 2 * exp(estParams(11));         % 消费测量误差
unoise_vol  = exp(estParams(11)); % 失业:误差波动率


fprintf('\n================ 基本参数结果 ================\n');
fprintf('Log-Likelihood: %.4f\n', LogL);
fprintf('----------------------------------------------\n');
fprintf('[长周期 ~7年] (Commercial Cycle)\n');
fprintf('Damping Factor (Long): %.4f  (越接近1越持久；约束限制)\n', rho_long);
fprintf('Cycle Period (Long):    %.2f Months (Target: 48-96)\n', period_long);
fprintf('Cycle Vol (Long):         %.5f (观察与趋势项大小对比)\n', scyc_long);
fprintf('Lambda (Long):             %.4f  (自由估计)\n', lambda_long);
fprintf('Okun Coef (Long):         %.4f (Expected < 0)\n', theta_long);
fprintf('----------------------------------------------\n');
fprintf('[短周期 ~2年] (Kitchin Cycle)\n');
fprintf('Damping Factor (Short): %.4f  (越接近1越持久；约束限制)\n', rho_short);
fprintf('Cycle Period (Short):     %.2f Months (Target: 12-32)\n', period_short);
fprintf('Cycle Vol (Short):          %.5f (观察与趋势项大小对比)\n', scyc_short);
fprintf('Lambda (Short):               %.4f  (自由估计)\n', lambda_short);
fprintf('Okun Coef (Short):          %.4f (Expected < 0)\n', theta_short);
fprintf('----------------------------------------------\n');

fprintf('Fiscal Sensitivity (Gamma):   %.4f\n', gamma_val);
fprintf('Fiscal unRate (Gamma):   %.4f\n', gamma_u_val);
fprintf('Income Noise (s_uy):               %.5f (锁死)\n', 0.00001);
fprintf('Consump Noise(s_uc):              %.5f (自由估计)\n', uc_val);
fprintf('UnempolyRate Noise (s_uu):    %.5f(自由估计)\n', unoise_vol);
fprintf('=================================================================\n');

fprintf('1. 假设检验：财政刺激显著性 (Fiscal Factor):\n');
% 提取 Gamma 的时变方差
gamma_cov_seq = squeeze(SmoothedStatesCov(5, 5, :)); 
gamma_se_final = sqrt(gamma_cov_seq(end)); % 取最终收敛的标准误
gamma_t_stat = gamma_val / gamma_se_final;
gamma_p_val = 2 * (1 - tcdf(abs(gamma_t_stat), T_len - 1)); % 自由度近似为 T
fprintf('Fiscal Sensitivity (Gamma):  %.4f\n', gamma_val);
fprintf('Fiscal Gamma SE:     %.4f\n', gamma_se_final);
fprintf('Fiscal Gamma T-Stat: %.4f\n', gamma_t_stat);
if gamma_p_val < 0.05
    fprintf('结论: 财政刺激效应显著 (Gamma != 0)\n');
else
    fprintf('结论: 财政刺激效应不显著\n');
end
fprintf('----------------------------------------------\n');
fprintf('2. 假设检验：Lambda_long > Lambda Short ?\n');
% 计算差值的标准误: Var(A-B) = Var(A) + Var(B) - 2*Cov(A,B)
idx_L = 4;  % lambda_long 在 params 中的索引
idx_S = 9; % lambda_short 在 params 中的索引
var_L = EstCov(idx_L, idx_L);
var_S = EstCov(idx_S, idx_S);
cov_LS = EstCov(idx_L, idx_S);

diff_lambda = lambda_long - lambda_short;
se_diff = sqrt(var_L + var_S - 2*cov_LS);
t_stat_diff = diff_lambda / se_diff;
p_val_diff = 1 - tcdf(t_stat_diff, T_len - 12); % 单侧检验 H0: L <= S
fprintf('  Diff (Long - Short): %.4f\n', diff_lambda);
fprintf('  T-Statistic:         %.4f\n', t_stat_diff);
fprintf('  P-Value (One-side):  %.4f\n', p_val_diff);
if p_val_diff < 0.05
    fprintf('  结论: 显著! 长周期收入被视为"持久收入"的程度更高。\n');
else
    fprintf('  结论: 不显著。长短周期敏感度无明显差异。\n');
end
fprintf('----------------------------------------------\n');
fprintf('\n================ Step 4 消费驱动力是什么？方差分解结果 ================\n');
% % --- 方差分解：消费驱动力分解 ---
% 1. 市场收入驱动部分
c_market_comp = lambda_long * psi_long_est + lambda_short * psi_short_est;
var_market = var(c_market_comp);
% 2. 财政转移驱动部分 (传入 Transfer_Factor)
% 假设 Transfer_Factor 在工作区中
c_fiscal_comp = gamma_val * trans_factor_adj; 
var_fiscal = var(c_fiscal_comp);
% 3. 忽略协方差项的简化对比
total_explained_var = var_market + var_fiscal;
fprintf('\n---- 消费波动驱动力对比 (Variance Contribution) ----\n');
fprintf('Market Income Contribution: %.2f%%\n', (var_market / total_explained_var)*100);
fprintf('Fiscal Transfer Contribution: %.2f%%\n', (var_fiscal / total_explained_var)*100);
fprintf('==================================================\n');

% =========================================================================
%% 6. 绘图
% =========================================================================

figure('Name', 'Step 4: Dual Cycle Results', 'Color', 'w');

% 图1: 两个周期的提取结果
subplot(3,2,1);
plot(data_time, psi_long_est, 'b-', 'LineWidth', 1.5); hold on;
plot(data_time, psi_short_est, 'color', [0.5 0.5 0.5], 'LineWidth', 1);
yline(0, 'k--');
title('Cycle Extraction: Long (Blue) vs. Short (Gray)');
legend(['Long Income Cycle (' num2str(period_long,'%.1f') 'm)'], ...
       ['Short Income Cycle (' num2str(period_short,'%.1f') 'm)']);
grid on; axis tight;

% 图2: 原始数据 vs 平滑趋势
subplot(3,2,2);
plot(data_time, ln_y_adj, 'k-', 'LineWidth', 2); hold on;
plot(data_time, ln_c_adj, 'b--', 'LineWidth', 2); 
plot(data_time, psi_long_est, 'r-', 'LineWidth', 1); 
title('Cycle Income vs. Cycle Consumption vs. Long Cycle');
legend('Income Cycle', 'Consumption Cycle', 'Long Cycle'); grid on; axis tight;

% 图3: 周期项与危机状态
subplot(3,2,3);
area(data_time, S_normalized * max(abs(u_gap_h)), 'FaceColor', [1 0.9 0.9], 'EdgeColor', 'none'); hold on;
plot(data_time, u_gap_h, 'k-', 'LineWidth', 1.5);
explained_u_gap = theta_long * psi_long_est + theta_short * psi_short_est;
plot(data_time, explained_u_gap, 'g--', 'LineWidth', 1.5);
yline(0, 'k--');
title('u_gap vs. explained Gap');
legend('Crisis peorid', 'Unemployment gap', 'Part expalined by Cycle');
grid on; axis tight;

% 图4: 隐含的周期传导
% c_cycle 和 lambda * y_cycle 的对比
c_cycle_implied_L = lambda_long .* psi_long_est; % 模型认为由收入传导过来的Long消费周期
c_cycle_implied_S = lambda_short .* psi_short_est;
c_cycle_fiscal = gamma_val * trans_factor_adj;
subplot(3,2,4);
% 这里粗略用 c - tau 作为消费周期观测值
plot(data_time, ln_c_adj, 'r-', 'LineWidth', 1); hold on;
plot(data_time, c_cycle_implied_L, 'b--', 'LineWidth', 1.5); 
plot(data_time, c_cycle_implied_S, 'g--', 'LineWidth', 1.5); 
plot(data_time, c_cycle_fiscal, 'k:', 'LineWidth', 1.5);
title('Consumption Cycle: Observed vs. Explained by Income Cycle');
legend('Total C Cycle (c - \tau)', 'Long-Income-Driven', 'Short-Income-Driven', 'Fiscal-Driven');
grid on; axis tight;

% 图5: 奥肯法则拟合
subplot(3,2,5);
plot(data_time, unRate, 'k-', 'LineWidth', 1); hold on;
plot(data_time, u_gap_h, 'r--', 'LineWidth', 1.5);
plot(data_time, u_trend_h, 'm--', 'LineWidth', 1.5);
title('Unemployment Fit');
legend('unRate', 'unRate (cycle)', 'Natural unRate');
grid on; axis tight;

% 图6: 消费周期到底解释了多少？
subplot(3,2,6);
% 这里粗略用 c - tau 作为消费周期观测值
plot(data_time, ln_c_adj, 'r-', 'LineWidth', 1); hold on;
plot(data_time, c_cycle_implied_L + c_cycle_implied_S + c_cycle_fiscal, 'b--', 'LineWidth', 1.5); 
title('Consumption Cycle: Observed vs. Explained by Income Cycle');
legend('Total C Cycle', 'Explained C Cycle');
grid on; axis tight;

% =========================================================================
%% 映射函数
% =========================================================================
function [A, B, C, D, Mean0, Cov0, StateType] = my_double_cycle_mapping(params, Transfer_Factor, T)
    
    % --- 参数解包 ---
    % 1-9 同原模型

    % 长周期 (Long)
    rho_L_min = 0.5; 
    rho_L_max = 0.995;
    rho_L = rho_L_min + (rho_L_max - rho_L_min) / (1 + exp(-params(1)));
    freq_L_min=2*pi/96;
    freq_L_max=2*pi/60;
    lam_L = freq_L_min + (freq_L_max - freq_L_min) / (1 + exp(-params(2)));
    scyc_L = exp(params(3));
    lambda_L = params(4);
    theta_L = -exp(params(5)); % Long term Okun Coefficient
    
    % 短周期 (Short)
    rho_S_min = 0.1; 
    rho_S_max = 0.95;
    rho_S = rho_S_min + (rho_S_max - rho_S_min) / (1 + exp(-params(6)));
    freq_S_min=2*pi/24;
    freq_S_max=2*pi/8;
    lam_S = freq_S_min + (freq_S_max - freq_S_min) / (1 + exp(-params(7)));
    scyc_S = exp(params(8));
    lambda_S = params(9);
    theta_S = params(10); % Short term Okun Coefficient

    % 噪声
    s_uc = exp(params(11)); % consumption noise
    %s_uu = exp(params(12)); % Unemplyment Noise

    
    % --- 2. 构造矩阵 (6 States) ---
    % States: [psi_L; psi*_L; psi_S; psi*_S; gamma]
    A = zeros(6,6);
    A(1,1) = rho_L*cos(lam_L); A(1,2) = rho_L*sin(lam_L);  % Long Cycle
    A(2,1) = -rho_L*sin(lam_L); A(2,2) = rho_L*cos(lam_L);
    A(3,3) = rho_S*cos(lam_S); A(3,4) = rho_S*sin(lam_S);  % Short Cycle
    A(4,3) = -rho_S*sin(lam_S); A(4,4) = rho_S*cos(lam_S);
    A(5,5) = 1; % gamma_c
    A(6,6) = 1; % gamma_u

    
    % B Matrix (6x7 Shocks)
    B = zeros(6,7);
    B(1,1) = scyc_L; % Long Real
    B(2,2) = scyc_L; % Long Imag
    B(3,3) = scyc_S; % Short Real
    B(4,4) = scyc_S; % Short Imag
    % Gamma 没有冲击 (Constant)

    % C Matrix (3 Obs x 6 States)
    C = zeros(3,6,T);
    C(1,1,:) = 1; C(1,3,:) = 1;                         % Obs 1: y = psi_L + psi_S (纯周期叠加)
    C(2,1,:) = lambda_L;C(2,3,:) = lambda_S; % Obs 2: c = lambda_L*psi_L + lambda_S*psi_S + gamma*Tr
    C(2,5,:) = reshape(Transfer_Factor, [1,1,T]); 
    C(3,1,:) = theta_L; C(3,3,:) = theta_S;    % Obs 3: u_gap = theta_L*psi_L + theta_S*psi_S
    C(3,6,:) = reshape(Transfer_Factor, [1,1,T]); 

    % D Matrix (3 Obs x 3 Obs Shocks)
    D = zeros(3,7);
    D(1,5) = 0.00001; % y_noise (Locked)
    D(2,6) = 2*s_uc;    % c_noise
    D(3,7) = s_uc;    % u_noise
    
    % Init
    Mean0 = zeros(6,1); Mean0(5) = 0.5; Mean0(6) = -0.5;
    Cov0 = diag([1, 1, 0.5, 0.5, 10, 10]);
    StateType = [0; 0; 0; 0; 2; 2];
end


% =========================================================================
% 辅助函数：计算负对数似然 (供 fminsearch 使用)
% =========================================================================
function nlogL = get_neg_logL(params, Data, Tr_vec, T)
    try
        [A, B, C, D, Mean0, Cov0, StateType] = my_double_cycle_mapping(params, Tr_vec, T);
        TempMdl = ssm(A, B, C, D, 'Mean0', Mean0, 'Cov0', Cov0, 'StateType', StateType);
        [~, logL] = filter(TempMdl, Data);
        nlogL = -logL;
        
        % 数值保护：防止 NaN 或 Inf 导致优化中断
        if isnan(nlogL) || isinf(nlogL)
            nlogL = 1e8; 
        end
    catch
        % 如果矩阵构造失败，返回大惩罚值，自动换方向，不报错停止
        nlogL = 1e8; 
    end
end
