
%% Export publication figures from the dual-cycle model workspace
project_root = fileparts(fileparts(mfilename('fullpath')));
%% ========================================================================
%  论文专用绘图脚本: 导出高清插图 (Publication Quality Plots)
%  前提: 请先运行完 Step 5 主程序，确保工作区里有 psi_long, data_time 等变量
% =========================================================================

% --- 1. 全局样式设置 (符合学术标准) ---
fig_width = 800;   % 图片宽度
fig_height = 500;  % 图片高度
font_name = 'Times New Roman'; % 论文标准字体
font_size = 12;    % 字号
line_width_thick = 1.5; % 主要线条粗细
line_width_thin = 1.0;  % 次要线条粗细

% 确保输出文件夹存在
%% ========================================================================
%  Paper Plot Generator: 批量导出 6 张独立高清图
%  说明：此脚本将原本的 subplot 拆解为独立的 figure，并按学术标准美化
% =========================================================================

% --- 0. 全局样式设置 (一次修改，全局生效) ---
save_dir = fullfile(project_root, 'results', 'figures');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

std_font = 'Times New Roman';     % 论文标准字体
std_size = 12;                    % 字号
line_w_main = 1.5;                % 主线宽
line_w_sub = 1.0;                 % 辅线宽
fig_pos = [100, 100, 600, 400];   % 图片默认尺寸 [x y w h]

fprintf('开始生成论文插图...\n');

%% --- 图 1: 周期提取结果 (Cycle Extraction) ---
f1 = figure('Color', 'w', 'Position', fig_pos, 'Name', 'Fig1_Cycles');
hold on;
yline(0, 'k-', 'HandleVisibility', 'off'); % 0轴
p_short = plot(data_time, psi_short_est, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', line_w_sub);
p_long  = plot(data_time, psi_long_est, 'b-', 'LineWidth', line_w_main);

title('Cycle Extraction: Long vs. Short', 'FontName', std_font, 'FontSize', std_size+1);
ylabel('Log Deviation', 'FontName', std_font);
% 注意：如果你已改为季度数据，建议将 'm' 改为 ' Qtrs'
legend([p_long, p_short], ...
       {['Long Cycle (' num2str(period_long,'%.1f') ')'], ...
        ['Short Cycle (' num2str(period_short,'%.1f') ')']}, ...
       'Location', 'best', 'Box', 'off', 'FontName', std_font);
axis tight; grid on; set(gca, 'FontName', std_font, 'FontSize', std_size);

% 导出
exportgraphics(f1, fullfile(save_dir, 'Fig1_Cycles.pdf'), 'ContentType', 'vector');
exportgraphics(f1, fullfile(save_dir, 'Fig1_Cycles.png'), 'Resolution', 300);


%% --- 图 2: 原始数据 vs 趋势 (Raw vs Trend) ---
f2 = figure('Color', 'w', 'Position', fig_pos, 'Name', 'Fig2_Trends');
hold on;
plot(data_time, ln_y_adj, 'k-', 'LineWidth', line_w_main);
plot(data_time, ln_c_adj, 'b--', 'LineWidth', line_w_main);
plot(data_time, psi_long_est, 'r-', 'LineWidth', line_w_sub);

title('Income vs. Consumption vs. Long Cycle', 'FontName', std_font, 'FontSize', std_size+1);
legend('Income (deTrended)', 'Consumption (deTrended)', 'Long Cycle', ...
       'Location', 'best', 'Box', 'off', 'FontName', std_font);
axis tight; grid on; set(gca, 'FontName', std_font, 'FontSize', std_size);

exportgraphics(f2, fullfile(save_dir, 'Fig2_Trends.pdf'), 'ContentType', 'vector');
exportgraphics(f2, fullfile(save_dir, 'Fig2_Trends.png'), 'Resolution', 300);


%% --- 图 3: 失业率缺口与危机 (Unemployment Gap) ---
f3 = figure('Color', 'w', 'Position', fig_pos, 'Name', 'Fig3_Unemployment');
hold on;
% 绘制阴影 (Crisis Area)
y_limits = [min(u_gap_h), max(u_gap_h)] * 1.1;
area_h = area(data_time, S_normalized * y_limits(2), 0, ...
    'FaceColor', [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
yline(0, 'k--', 'HandleVisibility', 'off');

% 绘制线条
p_ugap = plot(data_time, u_gap_h, 'k-', 'LineWidth', line_w_main);
% 计算解释部分
explained_u_gap = theta_long * psi_long_est + theta_short * psi_short_est; 
% 如果是 Unified Okun 模型，改为: theta_common * (psi_long_est + psi_short_est);
p_expl = plot(data_time, explained_u_gap, 'g--', 'LineWidth', line_w_main);

title('Unemployment Gap vs. Explained Gap', 'FontName', std_font, 'FontSize', std_size+1);
ylim(y_limits); % 限制Y轴
legend([area_h, p_ugap, p_expl], ...
       {'Crisis Period', 'Unemployment Gap', 'Explained by Cycles'}, ...
       'Location', 'best', 'Box', 'off', 'FontName', std_font);
axis tight; grid on; set(gca, 'FontName', std_font, 'FontSize', std_size);

exportgraphics(f3, fullfile(save_dir, 'Fig3_Unemployment.pdf'), 'ContentType', 'vector');
exportgraphics(f3, fullfile(save_dir, 'Fig3_Unemployment.png'), 'Resolution', 300);


%% --- 图 4: 消费周期驱动力分解 (Drivers of Consumption) ---
f4 = figure('Color', 'w', 'Position', fig_pos, 'Name', 'Fig4_Cons_Decomp');
hold on;
yline(0, 'k-', 'HandleVisibility', 'off');
% 计算分项
c_cycle_implied_L = lambda_long .* psi_long_est;
c_cycle_implied_S = lambda_short .* psi_short_est;
c_cycle_fiscal    = gamma_val * trans_factor_adj;

p_obs  = plot(data_time, ln_c_adj, 'r-', 'LineWidth', 1.2); 
p_long = plot(data_time, c_cycle_implied_L, 'b--', 'LineWidth', line_w_main); 
p_short= plot(data_time, c_cycle_implied_S, 'g--', 'LineWidth', line_w_sub); 
p_fisc = plot(data_time, c_cycle_fiscal, 'k:', 'LineWidth', line_w_main);

title('Decomposition of Consumption Cycle', 'FontName', std_font, 'FontSize', std_size+1);
legend([p_obs, p_long, p_short, p_fisc], ...
       {'Total C Cycle', 'Long-Income-Driven', 'Short-Income-Driven', 'Fiscal-Driven'}, ...
       'Location', 'best', 'Box', 'off', 'FontName', std_font);
axis tight; grid on; set(gca, 'FontName', std_font, 'FontSize', std_size);

exportgraphics(f4, fullfile(save_dir, 'Fig4_Cons_Decomp.pdf'), 'ContentType', 'vector');
exportgraphics(f4, fullfile(save_dir, 'Fig4_Cons_Decomp.png'), 'Resolution', 300);


%% --- 图 5: 奥肯法则拟合 (Okun Fit) ---
f5 = figure('Color', 'w', 'Position', fig_pos, 'Name', 'Fig5_Okun_Fit');
hold on;
plot(data_time, unRate, 'k-', 'LineWidth', 1.0); 
plot(data_time, u_gap_h, 'r--', 'LineWidth', line_w_main);
plot(data_time, u_trend_h, 'm--', 'LineWidth', line_w_main);

title('Unemployment Rate Decomposition', 'FontName', std_font, 'FontSize', std_size+1);
legend('Total UnRate', 'Cycle Component', 'Trend (by hamilton filter)', ...
       'Location', 'best', 'Box', 'off', 'FontName', std_font);
axis tight; grid on; set(gca, 'FontName', std_font, 'FontSize', std_size);

exportgraphics(f5, fullfile(save_dir, 'Fig5_Okun_Fit.pdf'), 'ContentType', 'vector');
exportgraphics(f5, fullfile(save_dir, 'Fig5_Okun_Fit.png'), 'Resolution', 300);


%% --- 图 6: 模型总解释力度 (Total Explained) ---
f6 = figure('Color', 'w', 'Position', fig_pos, 'Name', 'Fig6_Total_Explained');
hold on;
yline(0, 'k-', 'HandleVisibility', 'off');
% 总解释部分
total_explained = c_cycle_implied_L + c_cycle_implied_S + c_cycle_fiscal;

p_obs = plot(data_time, ln_c_adj, 'r-', 'LineWidth', line_w_sub); 
p_fit = plot(data_time, total_explained, 'b--', 'LineWidth', line_w_main); 

title('Consumption: Observed vs. Model Explained', 'FontName', std_font, 'FontSize', std_size+1);
legend([p_obs, p_fit], ...
       {'Observed Cycle', 'Model Explained (Income + Fiscal)'}, ...
       'Location', 'best', 'Box', 'off', 'FontName', std_font);
axis tight; grid on; set(gca, 'FontName', std_font, 'FontSize', std_size);

exportgraphics(f6, fullfile(save_dir, 'Fig6_Total_Explained.pdf'), 'ContentType', 'vector');
exportgraphics(f6, fullfile(save_dir, 'Fig6_Total_Explained.png'), 'Resolution', 300);

fprintf('所有图片导出完成！请查看文件夹: %s\n', save_dir);
