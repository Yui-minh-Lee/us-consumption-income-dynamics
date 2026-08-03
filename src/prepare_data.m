%% 第一阶段：数据准备与协整关系检查
clear; clc; close all;

project_root = fileparts(fileparts(mfilename('fullpath')));

% =========================================================================
% 1. 数据加载 (Load Data)
% =========================================================================
% 第一列是日期，第二列是 GDP(y)，第三列是 Consumption(c)
filename = fullfile(project_root, 'data', 'US_Data.csv');
data = readtable(filename); 

data = data(1:733, :); % 1959-01 to 2020-01; excludes the pandemic period.

% 取对数
ln_y = log(data.Real_income_excptTransfer);             % W875RX1, 实际市场收入
ln_c = log((data.Consumption ./ data.Price) * 100);  % 个人名义消费，转换实际变量
ln_dpi = log(data.Income ./ data.Price * 100);           % 个人名义可支配收入，转换实际值
ln_trans = log(data.Income ./ data.Price * 100 - data.Real_income_excptTransfer); % 实际可支配收入 - 实际市场收入 = 净转移支付
trans_factor = ln_dpi - ln_y;   % 转移支付因子

data_time = data.Date;

% 读完数据就好了，后面的代码不用管


%% =========================================================================
% 2. 协整直观检查：绘制 APC 路径
% ========================================================================
apc_ratio = ln_c - ln_y; % 对数消费-收入比 (log(C/Y))

figure('Name', 'Cointegration Check: c - y', 'Color', 'w');
plot(data_time, apc_ratio, 'LineWidth', 2, 'Color', 'b');
yline(mean(apc_ratio), '--r', 'Mean'); % 画出均值线
title('Log Consumption-Income Ratio (c_t - y_t)');
xlabel('Time'); ylabel('c - y');
grid on;
legend('c_t - y_t', 'Unconditional Mean');


%% =========================================================================
% 3. 构建 D_break
% ========================================================================
% 假设你在图中看到 2008Q4 发生了结构性断裂
T_len = length(ln_y);

% A. 结构突变 Dummy (Step Dummy)
break_date = datetime('2008-01-01');
 D_slope_2008 = zeros(T_len, 1);
idx_break = find(data_time >= break_date);
if ~isempty(idx_break)
     D_slope_2008(idx_break) = (1:length(idx_break))'; 
end

% B. 线性时间趋势捕捉 1975-2008 年爬升
Time_Trend = (1:T_len)'; 

% C. 斜率突变：1980 信贷扩张起点
% 允许 1980 之后的斜率发生改变
break_date_1980 = datetime('1980-01-01');
D_slope_1980 = zeros(T_len, 1);
idx_1980 = find(data_time >= break_date_1980);
if ~isempty(idx_1980)
    D_slope_1980(idx_1980) = (1:length(idx_1980))';
end



%% =========================================================================
% 4. 验证回归：残差是否平稳
% ========================================================================
% OLS 回归剔除结构突变的影响
X = [ones(T_len, 1), Time_Trend, D_slope_1980, D_slope_2008];

beta = (X'*X) \ (X' * apc_ratio);
residuals = apc_ratio - X * beta;

% 绘制剔除断点后的残差
figure('Name', 'Adjusted Residuals (With Trend)', 'Color', 'w');
subplot(2,1,1);
plot(data_time, apc_ratio, 'k--'); hold on;
plot(data_time, X*beta, 'r-', 'LineWidth', 2);
title('Model Fit: Including Linear Trend & COVID Dummy');
legend('Original Data', 'Fitted Step Break');
grid on;

subplot(2,1,2);
plot(data_time, residuals, 'b-', 'LineWidth', 1.5);
yline(0, 'k-');
title('Adjusted Residuals (Should be Stationary)');
grid on;

% 进行 ADF 单位根检验
% 原假设 H0: 存在单位根 (不平稳)
% 备择假设 H1: 平稳
[h, pValue, stat] = adftest(residuals, 'model', 'AR');

fprintf('======================================================\n');
fprintf('修正后的 ADF 检验结果:\n');
fprintf('p-value: %.4f\n', pValue);
if h == 1
    fprintf('结论: 拒绝原假设 (H0)。残差平稳！---> 现在的 D_break + Trend 设置有效。\n');
else
    fprintf('结论: 仍无法拒绝原假设 (H0)。\n');
    fprintf('建议: 可能需要 "Broken Trend" (趋势斜率发生改变)，而不仅仅是截距改变。\n');
end
fprintf('======================================================\n');

%% =========================================================================
% 5. 保存数据供后续使用
% =========================================================================
% 将处理好的 D_break 和标准化后的数据保存
processed_dir = fullfile(project_root, 'data', 'processed');
if ~exist(processed_dir, 'dir'), mkdir(processed_dir); end
save(fullfile(processed_dir, 'prepared_data.mat'), 'ln_y', 'ln_c', ...
    'trans_factor', 'Time_Trend', 'D_slope_1980', 'D_slope_2008', 'data');

X_trend = [ones(T_len,1), (1:T_len)'];
y_trend = X_trend * (((X_trend' * X_trend) \ X_trend' * ln_y));
% 2.清理 c 线性趋势、整理协整关系
apc_ratio = ln_c - ln_y;
X_det = [ones(T_len,1), Time_Trend, D_slope_1980, D_slope_2008];
% 计算协整关系的分段线性漂移趋势
coint_trend = X_det * ((X_det' * X_det) \ (X_det' * apc_ratio));
% 3.清理转移支付趋势
trans_factor_trend = X_trend * (((X_trend' * X_trend) \ X_trend' * trans_factor));
% 4.组装结果
ln_y_adj = ln_y -  y_trend;
ln_c_adj = ln_c - y_trend - coint_trend;
trans_factor_adj = trans_factor - trans_factor_trend;
