%% Step 5: 三变量时变参数模型
% 结构: 
%   y_t = tau_y + psi_t
%   c_t = tau_y + lambda_t * psi_t + gamma * Tr_t
%   u_gap_t = theta * psi_t  (直接观测失业率缺口)
%   fiscal_cycle
% 用失业率作为"锚"，强行纠正 2008 年的周期识别错误
% Plus:
%   1. 使用 "市场收入" (W875RX1) 识别纯粹的经济周期。
%   2. 将 "转移支付" (Trans_t) 独立出来解释消费，避免模型将财政刺激误判为趋势收入。
%   3. 现在y是市场收入 W875RX1。
%   4. 失业率趋势在模型外通过 HF 滤波移除，不再作为内部状态估计。

clear; clc; close all;

project_root = fileparts(fileparts(mfilename('fullpath')));

% =========================================================================
% 1. 数据准备
% =========================================================================
load(fullfile(project_root, 'data', 'processed', 'prepared_data.mat'));

% ---------- A. Hamilton滤波清洗：去趋势、整理协整关系 ----------
lead_len = 24;
lag_len = 12;
% 首部会有缺失, 截断起始NaN
data_time = data.Date(lead_len + lag_len : end);
T_len = length(data_time);

fprintf('正在使用 Hamilton 滤波进行去趋势处理 (LeadLength = %d, LagLength = %d)...\n', lead_len, lag_len);

% 1. 收入，不用去趋势 (Market Income)
ln_y_adj = ln_y(lead_len + lag_len : end);

% 2. 消费，去APC制度性漂移
apc_ratio = ln_c - ln_y;
[apc_trend_h, apc_cycle_h] = hfilter(apc_ratio, 'LeadLength', lead_len, 'LagLength', lag_len);
ln_c_adj = ln_c - apc_trend_h;
ln_c_adj = ln_c_adj(lead_len + lag_len : end);

% 3. 财政因子去趋势 (只保留刺激脉冲)
[trans_trend_h, trans_cycle_h] = hfilter(trans_factor, 'LeadLength', lead_len, 'LagLength', lag_len);
trans_factor_adj = trans_cycle_h(lead_len + lag_len : end);
fprintf('数据清洗完成。财政因子波动范围: %.4f ~ %.4f\n', min(trans_factor_adj), max(trans_factor_adj));


% 4. 失业率去趋势 HP滤波
unRate = data.Unemploy_rate / 100;
[u_trend_h, u_gap_h] = hfilter(unRate, 'LeadLength', lead_len, 'LagLength', lag_len);
unRate = unRate(lead_len + lag_len : end);
u_trend_h = u_trend_h(lead_len + lag_len : end);
u_gap_h = u_gap_h(lead_len + lag_len : end);

% 5.组装SSM输入数据
Data = [ln_y_adj, ln_c_adj, u_gap_h];


% ---------- B. 准备外部状态变量 S --------
% 这里的 S 必须是 "t-1" 时刻的信息，或者是当期外生的金融指标
S_raw = data.Unemploy_rate(lead_len + lag_len : end);
risk_threshold = mean(S_raw) + 1 * std(S_raw);  % threshold：1 std
S_normalized = double(S_raw > risk_threshold); % dummy array
S_normalized = [0; S_normalized(1:end-1)];       % lag 1

fprintf('状态变量 S_{t-1} 已构建并标准化 (Mean=%.2f, Std=%.2f)\n', mean(S_normalized), std(S_normalized));
fprintf('危机状态 S_{t-1} 占比: %.1f%%\n', mean(S_normalized)*100);


% =========================================================================
% 2. 继承 Step 2 的参数作为初始值
% =========================================================================
val_rho      = 0.9042 ;    % damping
val_period   = 80.0;     % cycle period (months)
val_eta      =  0.00127;   % trend vol
val_scyc     = 0.00956;   % cycle vol
val_rhocorr  = -0.3;   % correlation
val_lambda   = 0.5063;    % sensitivity (const part)
val_suc      = 0.0020;   % consump noise
val_suy      = 0.00001;   % income noise
val_lambda_beta = 0.5; % lambda sensitivity
val_okun      = -0.4;     % 奥肯系数 (预期为负)
val_suu       = 0.0040;      % 失业率测量误差

% --- 逆变换：带边界保护---- 
% 将实际经济参数转回无约束空间，作为 fminunc 的起点

% 防止 log(0) 或 log(negative)
clip = @(x, min_v, max_v) max(min_v + 1e-6, min(max_v - 1e-6, x));
inv_sigmoid = @(x, min_v, max_v) -log((max_v - min_v)/(x - min_v) - 1);

% 1. eta
eta_min = 1e-5;
eta_max = 0.02;
p1_0 =inv_sigmoid(clip(val_eta, eta_min, eta_max), eta_min, eta_max);

% 2. rho (damping)
rho_min = 0.50;
rho_max = 0.995;
p2_0 = inv_sigmoid(clip(val_rho, rho_min, rho_max), rho_min, rho_max);

% 3. freq
freq_min = 2*pi/96;
freq_max = 2*pi/12;
p3_0 = inv_sigmoid(clip(2*pi / val_period, freq_min, freq_max), freq_min, freq_max);

% 4. s_cyc (exp)
p4_0 = log(val_scyc);

% 5. rho_corr (-1, 1)
p5_0 = inv_sigmoid(clip(val_rhocorr, -1, 1), -1, 1);

% 6. lambda (const)
p6_0 = val_lambda;

% 7. s_uc (exp)
p7_0 = log(val_suc);

% 8. beta (sensitivity to S)
p8_0 = val_lambda_beta;

% 9. okun theta
p9_0 = log(-val_okun);

% 10. U Noise (Log): > 0 constraint
p10_0 = log(val_suu);


params0 = [p1_0, p2_0, p3_0, p4_0, p5_0, p6_0, p7_0, p8_0, p9_0, p10_0];


% =========================================================================
% 3. 定义 TVP 模型
% =======================================================================
% 使用新的三变量映射函数
Mdl_Tri = ssm(@(p) my_trivariate_mapping(p, S_normalized, trans_factor_adj, T_len));


% =========================================================================
% 4. 模型估计 (改进版：fminsearch + fminunc 混合策略)
% =======================================================================
fprintf('开始 Step 4 双时变参数估计 (TVP-SSM)...\n');
% --- 阶段 A: 使用 fminsearch 单纯形法搜索 ---
fprintf('\n--- Phase 1: Pre-optimization with fminsearch (Simplex) ---\n');
ObjFun = @(p) get_neg_logL(p, Data, S_normalized, trans_factor_adj, T_len);  % 定义目标函数句柄 (负对数似然)
options_search = optimset('MaxFunEvals',50000,'MaxIter',5000);  % 设置 fminsearch 选项
t_start = tic; % 记录时间
[params_search, fval_search] = fminsearch(ObjFun, params0, options_search); % 运行 fminsearch
fprintf('fminsearch 完成! 耗时: %.1f 秒, Min NLogL: %.4f\n', toc(t_start), fval_search);

% --- 阶段 B: 将 fminsearch 的结果作为 fminunc 新起点 ---
fprintf('\n--- Phase 2: Polishing with estimate (fminunc) ---\n');
fprintf('使用 fminsearch 找到的最优解作为新起点，开始计算 Hessian 和标准误...\n');

params0 = params_search; % 更新初始值
options = optimoptions('fminunc', 'MaxFunctionEvaluations', 50000, ...
    'MaxIterations', 5000, 'Display', 'iter', ...
    'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6);

try
    [EstMdl, estParams, EstCov, logL] = ...
        estimate(Mdl_Tri, Data, params0, 'Options', options);
    fprintf('Step 4 最终估计完成！LogL: %.4f\n', logL);
catch ME
    fprintf('估计失败: %s\n', ME.message);
    return;
end

% =========================================================================
%% 5. 结果解析与绘图
% =======================================================================
% 提取状态
[SmoothedStates, LogL, Output] = smooth(EstMdl, Data);
SmoothedStatesCov = Output.SmoothedStatesCov;

tau_y_est = SmoothedStates(:, 1); % income trend
psi_est    = SmoothedStates(:, 2); % common cycle
gamma_est_seq = SmoothedStates(:, 4); % 常数序列

% 提取参数 lambda
lambda_const_est = estParams(6);
beta_est               = estParams(8);
beta_est_se          = sqrt(EstCov(8,8));
beta_t_stat   = beta_est / beta_est_se;                                 % 时变系数 beta ≠ 0 双侧检验
beta_p_val    = 2 * (1 - tcdf(abs(beta_t_stat), T_len - 11));
lambda_path = lambda_const_est + beta_est * S_normalized;  % 构造 lambda 时变路径

gamma_val  = mean(gamma_est_seq); % 财政乘数/敏感度

theta_est    = -exp(estParams(9)); % 奥肯系数
unoise_vol  = exp(estParams(10)); % 失业:误差波动率
scyc_val      = exp(estParams(4));  % 周期方差

% 1. 趋势波动 (映射回 [1e-5, 0.0005])
sig_eta_val = eta_min + (eta_max - eta_min) / (1 + exp(-estParams(1)));
% 2. 阻尼因子 (映射回 [0.90, 0.995])
rho_damp = rho_min + (rho_max - rho_min) / (1 + exp(-estParams(2)));
% 频率与周期长度
freq_val = freq_min + (freq_max - freq_min) / (1 + exp(-estParams(3)));
period_months = 2 * pi / freq_val; % 周期长度 (月)
% 4. 相关系数 rho_corr
rho_corr_val = 2 / (1 + exp(-estParams(5))) - 1; 
% 6. 消费测量误差
uc_val = exp(estParams(7));

fprintf('\n================ 基本参数结果 ================\n');
fprintf('Damping Factor:       %.4f  (越接近1越持久；约束限制)\n', rho_damp);
fprintf('Cycle Period:         %.2f Months (Target: 18-96)\n', period_months);
fprintf('Trend Vol (sig_eta):  %.5f (约束限制)\n', sig_eta_val);
fprintf('Cycle Vol (sig_cyc):  %.5f (观察与趋势项大小对比)\n', scyc_val);
fprintf('Correlation (rho):    %.4f  (自由估计)\n', rho_corr_val);
fprintf('Sensitivity(lambda): %.4f  (自由估计)\n', lambda_const_est);
fprintf('Fiscal Sensitivity (Gamma):  %.4f\n', gamma_val);
fprintf('Income Noise (s_uy):  %.5f (锁死)\n', 0.00001);
fprintf('Consump Noise(s_uc):  %.5f (自由估计)\n', uc_val);
fprintf('UnempolyRate Noise (s_uu):   %.5f(自由估计)\n', unoise_vol);
fprintf('=================================================================\n');

fprintf('\n================ Step 4 TVP 结果 ================\n');
fprintf('1. 敏感性渠道 (Sensitivity Channel):\n');
fprintf('Market Income Sensitivity: %.4f\n', lambda_const_est);
fprintf('Fiscal Sensitivity (Gamma):  %.4f\n', gamma_val);
fprintf('Beta (Coef of S):  %.4f\n', beta_est);
fprintf('Standard Error: %.4f\n', beta_est_se);
fprintf('T-Statistic:    %.4f\n', beta_t_stat);
fprintf('P-Value:        %.4f\n', beta_p_val);
if beta_p_val < 0.05
    fprintf('结论: 拒绝 H0 (***)。存在显著的时变过度敏感性！\n');
else
    fprintf('结论: 无法拒绝 H0。未发现显著的时变特征。\n');
end
fprintf('----------------------------------------------\n');
fprintf('2. 失业率相关 (UnRate):\n');
fprintf('Okun Coef:   %.4f (Expected < 0)\n', theta_est);
fprintf('----------------------------------------------\n');
fprintf('3. 财政刺激相关 (Fiscal Factor):\n');
% 提取 Gamma 的时变方差
gamma_cov_seq = squeeze(SmoothedStatesCov(4, 4, :)); 
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
fprintf('\n================ Step 4 消费驱动力是什么？方差分解结果 ================\n');
% % --- 方差分解：消费驱动力分解 ---
% 1. 市场收入驱动部分
c_market_comp = lambda_path .* psi_est;
var_market = var(c_market_comp);

% 2. 财政转移驱动部分 (注意要传入 Transfer_Factor)
% 假设 Transfer_Factor 在工作区中
c_fiscal_comp = gamma_val * trans_factor_adj; 
var_fiscal = var(c_fiscal_comp);

total_explained_var = var_market + var_fiscal; % 忽略协方差项的简化对比
fprintf('\n---- 消费波动驱动力对比 (Variance Contribution) ----\n');
fprintf('Market Income Contribution: %.2f%%\n', (var_market / total_explained_var)*100);
fprintf('Fiscal Transfer Contribution: %.2f%%\n', (var_fiscal / total_explained_var)*100);
fprintf('==================================================\n');

% =========================================================================
%% 6. 绘图
% =======================================================================

figure('Name', 'Step 4: TVP Results', 'Color', 'w');

% 图1: lambda 路径 + 外部状态变量 S
subplot(3,2,1);
plot(data_time, lambda_path, 'r-', 'LineWidth', 2); hold on;
plot(data_time, S_normalized, 'b-', 'LineWidth', 1); 
yline(0, 'k--');
title(['Time-Varying Sensitivity \lambda_t = \lambda_0 + \beta S_{t-1} (Beta = ' num2str(beta_est, '%.4f') ') plus S_t']);
ylabel('\lambda_t'); grid on; axis tight;

% 图2: 原始数据 vs 平滑趋势
subplot(3,2,2);
plot(data_time, ln_y_adj, 'k-', 'LineWidth', 2); hold on;
plot(data_time, ln_c_adj, 'b--', 'LineWidth', 2); 
plot(data_time, tau_y_est, 'r-', 'LineWidth', 1); 
title('Log Income (Adjusted) vs. log Consumption vs. Smooth Trend (\tau_t)');
legend('Income', 'Consumption', 'Smooth Trend'); grid on; axis tight;

% 图3: 周期项与危机状态
subplot(3,2,3);
area(data_time, S_normalized * max(abs(u_gap_h)), 'FaceColor', [1 0.9 0.9], 'EdgeColor', 'none'); hold on;
plot(data_time, u_gap_h, 'k-', 'LineWidth', 1.5);
explained_u_gap = theta_est * psi_est;
plot(data_time, explained_u_gap, 'r--', 'LineWidth', 1.5);
yline(0, 'k--');
title('u_gap vs. explained Gap');
legend('Crisis peorid', 'Unemployment gap', 'Part expalined by Cycle', 'Part explained by fiscal');
grid on; axis tight;

% 图4: 隐含的周期传导
% 展示 c_cycle 和 lambda * y_cycle 的对比
c_cycle_implied = lambda_path .* psi_est; % 模型认为由收入传导过来的消费周期
c_cycle_fiscal = gamma_val * trans_factor_adj;
subplot(3,2,4);
% 这里粗略用 c - tau 作为消费周期观测值
c_cycle_obs = ln_c_adj - tau_y_est; 
plot(data_time, c_cycle_obs, 'r-', 'LineWidth', 1); hold on;
plot(data_time, c_cycle_implied, 'b--', 'LineWidth', 1.5); 
plot(data_time, c_cycle_fiscal, 'k:', 'LineWidth', 1.5);
title('Consumption Cycle: Observed vs. Explained by Income Cycle');
legend('Total C Cycle (c - \tau)', 'Income-Driven (\lambda_t \psi_t)', 'Fiscal-Driven');
grid on; axis tight;

% 图5: 奥肯法则拟合
subplot(3,2,5);
plot(data_time, unRate, 'k-', 'LineWidth', 1); hold on;
plot(data_time, u_gap_h, 'r--', 'LineWidth', 1.5);
plot(data_time, u_trend_h, 'm--', 'LineWidth', 1.5);
title(['Unemployment Fit (Okun \theta = ' num2str(theta_est, '%.2f') ')']);
legend('unRate', 'unRate (cycle)', 'Natural unRate');
grid on; axis tight;


% =========================================================================
%% 映射函数
% =========================================================================
function [A, B, C, D, Mean0, Cov0, StateType] = my_trivariate_mapping(params, S_vec, Transfer_Factor, T)
    
    % --- 参数解包 ---
    eta_min=1e-5; 
    eta_max=0.02;
    s_eta = eta_min + (eta_max - eta_min) / (1 + exp(-params(1)));
    
    rho_min=0.5; 
    rho_max=0.995;
    rho = rho_min + (rho_max - rho_min) / (1 + exp(-params(2)));
    
    freq_min=2*pi/96; 
    freq_max=2*pi/12;
    lam_c = freq_min + (freq_max - freq_min) / (1 + exp(-params(3)));
    
    scyc = exp(params(4));
    rho_corr = 2 / (1 + exp(-params(5))) - 1;
    
    lam_const = params(6);
    s_uc = exp(params(7));
    beta = params(8);
    theta = -exp(params(9));
    s_uu = exp(params(10));
    
    % --- 构造向量 ---
    lambda_vec = reshape(lam_const + beta * S_vec, [1,1,T]);
    
    % --- 构造矩阵 (4 States) ---
    % States: [tau_y; psi; psi*; tau_u; fiscal]
    % A Matrix (4x4)
    A = zeros(4,4);
    A(1,1) = 1; % tau_y
    A(2,2) = rho*cos(lam_c); A(2,3) = rho*sin(lam_c); % cycle
    A(3,2) = -rho*sin(lam_c); A(3,3) = rho*cos(lam_c);
    A(4,4) = 1; % gamma_c
    
    % B Matrix (4x7 Shocks)
    B = zeros(4,7);
    corr_term = sqrt(max(0, 1 - rho_corr^2));
    B(1,1) = s_eta;                       % Shock 1: Trend_y
    B(2,1) = rho_corr * scyc;       
    B(2,2) = scyc * corr_term;     % Shock 2: Cycle (Corr)
    B(3,3) = scyc;                         % Shock 3: Cycle (Indep)
    % B(4,:) are all zeros -> Gamma is constant over time
    
    % C Matrix (3 Obs x 5 States)
    C = zeros(3,4,T);
    C(1,1,:) = 1; C(1,2,:) = 1;                 % Obs 1: y_t = tau_y + psi
    C(2,1,:) = 1; C(2,2,:) = lambda_vec; % Obs 2: c_t = tau_y + lambda_t * psi
    C(2,4,:) = reshape(Transfer_Factor, [1,1,T]); % Tr_t multiplies gamma
    C(3,2,:) = theta;                                % Obs 3: u_gap = theta * psi
   

    % D Matrix (3 Obs x 3 Obs Shocks)
    D = zeros(3,7);
    D(1,5) = 0.00001; % y_noise (Locked)
    D(2,6) = s_uc;    % c_noise
    D(3,7) = s_uu;    % u_noise
    
    % Init
    Mean0 = zeros(4,1);
    Mean0(4) = 0.5;
    Cov0 = diag([10, 1, 1, 10]); % Diffuse for trends
    StateType = [2; 0; 0; 2];
end


% =========================================================================
% 辅助函数：计算负对数似然
% =========================================================================
function nlogL = get_neg_logL(params, Data, S_vec, Tr_vec, T)
    try
        [A, B, C, D, Mean0, Cov0, StateType] = my_trivariate_mapping(params, S_vec, Tr_vec, T);
        TempMdl = ssm(A, B, C, D, 'Mean0', Mean0, 'Cov0', Cov0, 'StateType', StateType);
        [~, logL] = filter(TempMdl, Data);
        nlogL = -logL;
        
        % 数值保护：防止 NaN 或 Inf 导致优化中断
        if isnan(nlogL) || isinf(nlogL)
            nlogL = 1e8; 
        end
    catch
        % 如果矩阵构造失败（如非正定），返回大惩罚值，自动换方向，不报错停止
        nlogL = 1e8; 
    end
end
