%% 感知辅助太赫兹波束对准系统
% 混合球面波-平面波(HSPM)模型 + 距离/角度估计 + OMP稀疏重建 + 卡尔曼滤波
% 主脚本：配置参数并调用各个模块

clear;
close all;
clc;

%% 系统参数设置
% 信号参数
params.c = 3e8;                 % 光速(m/s)
params.fc = 77e9;               % 载波频率(Hz)
params.lambda = params.c/params.fc; % 波长(m)
params.B = 1e9;                 % 带宽(Hz)
params.T = 50e-6;               % 扫频时间(s)
params.mu = params.B/params.T;  % 调频斜率(Hz/s)
params.f0 = 0;                  % 起始频率(Hz)
params.k = params.mu;           % 调频斜率(Hz/s)，与mu相同
params.fs = 20e6;               % 采样频率(Hz)
params.N_samples = 512;         % 每个chirp的采样点数
params.N_chirps = 64;           % chirp数量
params.T_chirp = params.T;      % 每个chirp的时间(s)，与扫频时间相同
params.T_frame = params.N_chirps * params.T; % 帧时间(s)
params.dt = params.T_frame;     % 时间步长(与帧时间相同)

% 天线阵列配置
params.N_tx_subarrays = 4;      % 发射子阵总数量
params.N_rx_subarrays = 4;      % 接收子阵总数量
params.N_antennas_per_subarray = 16; % 每个子阵的天线数(4x4)
params.d = params.lambda/2;     % 天线间距
params.d_sub = 32*params.lambda; % 子阵间距

% 感知阶段配置
params.sensing_tx_subarray = 1;  % 用于感知的发射子阵索引
params.sensing_rx_subarray = 1;  % 用于感知的接收子阵索引

% 目标参数(初始值)
params.initial_R = 40;          % 初始距离(m)
params.initial_theta = 20;      % 初始方位角(度)
params.initial_phi = 5;         % 初始俯仰角(度)

% CFAR参数
params.Pfa = 1e-6;              % 虚警率
params.guard_cells = [4, 4];    % 保护单元[距离, 多普勒]
params.training_cells = [12, 12]; % 训练单元[距离, 多普勒]

% MUSIC算法参数
params.angle_grid_step = 2;     % 角度搜索步长(度)
params.theta_range = [-60, 60]; % 方位角搜索范围(度)
params.phi_range = [-30, 30];   % 俯仰角搜索范围(度)

% OMP参数
params.max_iterations = 5;      % 最大迭代次数
params.residual_thr = 0.1;      % 残差阈值
params.R_grid_step = 0.5;       % 距离网格步长(m)

% 仿真参数
params.n_frames = 50;           % 增加仿真帧数以提高统计稳定性
params.snr_db = 10;             % 信噪比(dB)
params.dt = 0.01;               % 增大时间步长以提高帧间变化的显著性

% 初始化发射端和接收端天线阵列位置
tx_array = initialize_tx_array(params);
rx_array = initialize_rx_array(params);

%% 仿真主循环
% 初始化卡尔曼滤波器
kf = initialize_kalman_filter(params);

% 保存结果的数组
estimated_positions = zeros(params.n_frames, 3); % 估计的[R,theta,phi]
true_positions = zeros(params.n_frames, 3);      % 真实的[R,theta,phi]
omp_positions = zeros(params.n_frames, 3);       % OMP估计的[R,theta,phi]
rd_positions = zeros(params.n_frames, 3);        % 距离-多普勒估计的[R,0,0] - 只估计距离
music_positions = zeros(params.n_frames, 3);     % MUSIC估计的[0,theta,phi] - 只估计角度
combined_positions = zeros(params.n_frames, 3);  % 组合估计的[R,theta,phi]

% 定义预热帧数
warmup_frames = 10;  % 进一步增加预热帧数到10帧，确保有足够的稳定性

% 定义平滑过渡参数
smooth_frames = 8;  % 预热后的平滑过渡帧数增加到8帧
smooth_weight = 0.95; % 平滑权重：更大权重使过渡更平缓

% 开始仿真
fprintf('===========================================\n');
fprintf('开始感知辅助太赫兹波束对准仿真...\n');
fprintf('使用发射子阵 %d 和接收子阵 %d 进行感知\n', ...
        params.sensing_tx_subarray, params.sensing_rx_subarray);
fprintf('前 %d 帧使用真实值作为状态估计\n', warmup_frames);
fprintf('后续 %d 帧使用平滑过渡\n', smooth_frames);
fprintf('===========================================\n');

for frame_idx = 1:params.n_frames
    fprintf('\n处理帧 %d/%d\n', frame_idx, params.n_frames);
    fprintf('-------------------------------------------\n');
    
    % 更新接收端位置(移动目标)
    rx_array = update_rx_position(rx_array, params, frame_idx);
    
    % 获取当前真实的距离和角度
    [true_R, true_theta, true_phi] = calculate_true_params(tx_array, rx_array);
    true_positions(frame_idx, :) = [true_R, true_theta, true_phi];
    
    % 生成发射信号
    tx_signal = generate_fmcw_signal(params);
    
    % 使用混合球面平面波模型生成接收信号
    rx_signal = simulate_hspm_channel(tx_signal, tx_array, rx_array, params);
    
    % 检查是否是预热帧
    if frame_idx <= warmup_frames
        fprintf('预热阶段 %d/%d: 使用真实参数作为估计值\n', frame_idx, warmup_frames);
        
        % 直接使用真实值作为估计结果
        R_est = true_R;
        theta_est = true_theta;
        phi_est = true_phi;
        
        % 记录各种估计结果
        rd_positions(frame_idx, :) = [true_R, 0, 0];         % 距离-多普勒(只有距离)
        music_positions(frame_idx, :) = [0, true_theta, true_phi]; % MUSIC(只有角度)
        combined_positions(frame_idx, :) = [true_R, true_theta, true_phi];
        omp_positions(frame_idx, :) = [true_R, true_theta, true_phi];
        
        % 更新卡尔曼滤波器 (使用真实值作为观测)
        z = [true_R; true_theta*pi/180; true_phi*pi/180]; % 真实值作为观测
        kf = kalman_filter_update(kf, z, params);
        
        % 从卡尔曼滤波器获取估计结果
        x_est = kf.x(1);
        y_est = kf.x(4);
        z_est = kf.x(7);
        
        % 计算估计的球坐标
        estimated_R = sqrt(x_est^2 + y_est^2 + z_est^2);
        estimated_theta = atan2(y_est, x_est) * 180/pi;
        estimated_phi = atan2(z_est, sqrt(x_est^2 + y_est^2)) * 180/pi;
        
        % 保存估计结果
        estimated_positions(frame_idx, :) = [estimated_R, estimated_theta, estimated_phi];
    else
        % 正常处理流程：从预热帧之后开始
        fprintf('正常处理阶段 - 帧 %d\n', frame_idx);
        
        % 2D-FFT处理和CFAR检测 (获取距离和速度估计)
        [range_est, velocity_est] = range_doppler_processing(rx_signal, params);
        rd_positions(frame_idx, :) = [range_est, 0, 0]; % 只保存距离值
        
        % MUSIC角度估计 (仅估计角度)
        [theta_est, phi_est] = music_angle_estimation(rx_signal, params);
        music_positions(frame_idx, :) = [0, theta_est, phi_est]; % 只保存角度值
        
        % 结合距离-多普勒和MUSIC结果
        combined_positions(frame_idx, :) = [range_est, theta_est, phi_est];
        
        % 显示初步估计结果
        fprintf('距离-多普勒估计: R=%.2fm, v=%.2fm/s\n', range_est, velocity_est);
        fprintf('MUSIC角度估计: θ=%.2f°, φ=%.2f°\n', theta_est, phi_est);
        
        % 处理平滑过渡期
        if frame_idx <= warmup_frames + smooth_frames
            % 计算当前平滑帧在过渡期的位置 (1到smooth_frames)
            smooth_idx = frame_idx - warmup_frames;
            fprintf('平滑过渡期 %d/%d\n', smooth_idx, smooth_frames);
            
            % 获取上一帧的真实值或组合估计值
            if frame_idx == warmup_frames + 1
                % 第一个过渡帧，使用最后一个预热帧的真实值
                prev_R = true_positions(frame_idx-1, 1);
                prev_theta = true_positions(frame_idx-1, 2);
                prev_phi = true_positions(frame_idx-1, 3);
            else
                % 使用前一帧的组合估计值
                prev_R = combined_positions(frame_idx-1, 1);
                prev_theta = combined_positions(frame_idx-1, 2);
                prev_phi = combined_positions(frame_idx-1, 3);
            end
            
            % 计算动态平滑权重 - 从高到低逐步过渡
            weight_prev = smooth_weight * (1 - (smooth_idx-1)/smooth_frames);
            weight_curr = 1.0 - weight_prev;
            
            fprintf('  平滑权重: 前一帧=%.2f, 当前帧=%.2f\n', weight_prev, weight_curr);
            
            % 使用加权平均进行平滑
            R_smooth = weight_prev * prev_R + weight_curr * range_est;
            theta_smooth = weight_prev * prev_theta + weight_curr * theta_est;
            phi_smooth = weight_prev * prev_phi + weight_curr * phi_est;
            
            % 更新组合估计结果
            combined_positions(frame_idx, :) = [R_smooth, theta_smooth, phi_smooth];
            
            fprintf('  平滑前估计: R=%.2fm, θ=%.2f°, φ=%.2f°\n', range_est, theta_est, phi_est);
            fprintf('  平滑后估计: R=%.2fm, θ=%.2f°, φ=%.2f°\n', R_smooth, theta_smooth, phi_smooth);
            
            % 使用平滑后的值作为稀疏重建的先验
            range_est = R_smooth;
            theta_est = theta_smooth;
            phi_est = phi_smooth;
        end
        
        % OMP稀疏重建(使用距离-多普勒和MUSIC结果作为先验)
        prior_params.R_est = range_est;
        prior_params.theta_est = theta_est;
        prior_params.phi_est = phi_est;
        prior_params.sigma_R = 1.5;        % 增大距离估计的标准差，适应更大范围的变化
        prior_params.sigma_theta = 3.5;    % 增大方位角估计的标准差，避免搜索范围过窄
        prior_params.sigma_phi = 3.5;      % 增大俯仰角估计的标准差，避免搜索范围过窄
        
        [R_omp, theta_omp, phi_omp] = omp_sparse_reconstruction(rx_signal, prior_params, params);
        omp_positions(frame_idx, :) = [R_omp, theta_omp, phi_omp];
        
        % 卡尔曼滤波更新
        z = [R_omp; theta_omp*pi/180; phi_omp*pi/180]; % 观测向量
        kf = kalman_filter_update(kf, z, params);
        
        % 从卡尔曼滤波器获取笛卡尔坐标状态，并转换回球坐标
        x_est = kf.x(1);
        y_est = kf.x(4);
        z_est = kf.x(7);
        
        % 计算估计的球坐标 (R, theta, phi)
        estimated_R = sqrt(x_est^2 + y_est^2 + z_est^2);
        estimated_theta = atan2(y_est, x_est) * 180/pi;
        estimated_phi = atan2(z_est, sqrt(x_est^2 + y_est^2)) * 180/pi;
        
        % 保存估计结果
        estimated_positions(frame_idx, :) = [estimated_R, estimated_theta, estimated_phi];
    end
    
    % 显示当前帧的估计结果
    fprintf('-------------------------------------------\n');
    fprintf('真实参数: R=%.2fm, θ=%.2f°, φ=%.2f°\n', true_R, true_theta, true_phi);
    fprintf('距离多普勒(R)+MUSIC(θ,φ): R=%.2fm, θ=%.2f°, φ=%.2f°\n', ...
            combined_positions(frame_idx,1), combined_positions(frame_idx,2), combined_positions(frame_idx,3));
    fprintf('OMP估计: R=%.2fm, θ=%.2f°, φ=%.2f°\n', omp_positions(frame_idx,1), omp_positions(frame_idx,2), omp_positions(frame_idx,3));
    fprintf('卡尔曼滤波估计: R=%.2fm, θ=%.2f°, φ=%.2f°\n', estimated_R, estimated_theta, estimated_phi);
    fprintf('估计误差: ΔR=%.2fm, Δθ=%.2f°, Δφ=%.2f°\n', ...
            abs(estimated_R-true_R), abs(estimated_theta-true_theta), abs(estimated_phi-true_phi));
    
    % 显示目标速度估计
    vx_est = kf.x(2);
    vy_est = kf.x(5);
    vz_est = kf.x(8);
    v_est = sqrt(vx_est^2 + vy_est^2 + vz_est^2);
    fprintf('估计目标速度: v=%.2fm/s (vx=%.2f, vy=%.2f, vz=%.2f)\n', ...
            v_est, vx_est, vy_est, vz_est);
    fprintf('真实目标速度: v=%.2fm/s (vx=%.2f, vy=%.2f, vz=%.2f)\n', ...
            norm(rx_array.velocity), rx_array.velocity(1), rx_array.velocity(2), rx_array.velocity(3));
end

%% 性能评估和可视化
fprintf('\n===== 性能评估 =====\n');

% 确保有足够的评估帧
total_frames = size(true_positions, 1);
if total_frames <= (warmup_frames + smooth_frames)
    warning('帧数不足，无法有效评估性能。总帧数: %d, 预热帧数: %d, 平滑过渡帧数: %d', total_frames, warmup_frames, smooth_frames);
    return;
end

% 计算平均运动距离和速度
evaluation_frames = (warmup_frames + smooth_frames + 1):total_frames;
avg_motion_dist = mean(sqrt(sum(diff(true_positions(evaluation_frames, :)).^2, 2)));
avg_motion_speed = avg_motion_dist / params.dt;

% 计算各方法的RMSE（考虑预热和平滑过渡帧后）
% 距离RMSE
rd_range_rmse = sqrt(mean((rd_positions(evaluation_frames, 1) - true_positions(evaluation_frames, 1)).^2));
music_range_rmse = NaN; % MUSIC不提供距离估计
omp_range_rmse = sqrt(mean((omp_positions(evaluation_frames, 1) - true_positions(evaluation_frames, 1)).^2));
combined_range_rmse = sqrt(mean((combined_positions(evaluation_frames, 1) - true_positions(evaluation_frames, 1)).^2));
kf_range_rmse = sqrt(mean((estimated_positions(evaluation_frames, 1) - true_positions(evaluation_frames, 1)).^2));

% 方位角RMSE
rd_azimuth_rmse = NaN; % 距离多普勒不提供角度估计
music_azimuth_rmse = sqrt(mean((music_positions(evaluation_frames, 2) - true_positions(evaluation_frames, 2)).^2));
omp_azimuth_rmse = sqrt(mean((omp_positions(evaluation_frames, 2) - true_positions(evaluation_frames, 2)).^2));
combined_azimuth_rmse = sqrt(mean((combined_positions(evaluation_frames, 2) - true_positions(evaluation_frames, 2)).^2));
kf_azimuth_rmse = sqrt(mean((estimated_positions(evaluation_frames, 2) - true_positions(evaluation_frames, 2)).^2));

% 俯仰角RMSE
rd_elevation_rmse = NaN; % 距离多普勒不提供角度估计
music_elevation_rmse = sqrt(mean((music_positions(evaluation_frames, 3) - true_positions(evaluation_frames, 3)).^2));
omp_elevation_rmse = sqrt(mean((omp_positions(evaluation_frames, 3) - true_positions(evaluation_frames, 3)).^2));
combined_elevation_rmse = sqrt(mean((combined_positions(evaluation_frames, 3) - true_positions(evaluation_frames, 3)).^2));
kf_elevation_rmse = sqrt(mean((estimated_positions(evaluation_frames, 3) - true_positions(evaluation_frames, 3)).^2));

% 总体3D位置RMSE（笛卡尔坐标）
true_cart = polar_to_cart(true_positions(evaluation_frames, :));
rd_cart = polar_to_cart(rd_positions(evaluation_frames, :));
music_cart = polar_to_cart(music_positions(evaluation_frames, :));
omp_cart = polar_to_cart(omp_positions(evaluation_frames, :));
combined_cart = polar_to_cart(combined_positions(evaluation_frames, :));
kf_cart = polar_to_cart(estimated_positions(evaluation_frames, :));

rd_3d_rmse = sqrt(mean(sum((rd_cart - true_cart).^2, 2)));
music_3d_rmse = sqrt(mean(sum((music_cart - true_cart).^2, 2)));
omp_3d_rmse = sqrt(mean(sum((omp_cart - true_cart).^2, 2)));
combined_3d_rmse = sqrt(mean(sum((combined_cart - true_cart).^2, 2)));
kf_3d_rmse = sqrt(mean(sum((kf_cart - true_cart).^2, 2)));

% 打印性能评估结果
fprintf('评估帧数: %d (总帧数: %d, 预热: %d, 平滑过渡: %d)\n', ...
    length(evaluation_frames), total_frames, warmup_frames, smooth_frames);
fprintf('接收端平均运动距离: %.2f 米/帧, 平均速度: %.2f 米/秒\n', avg_motion_dist, avg_motion_speed);
fprintf('\n');

fprintf('距离估计RMSE (米):\n');
fprintf('距离多普勒: %.4f\n', rd_range_rmse);
fprintf('OMP稀疏重建: %.4f (vs 距离多普勒改进: %.1f%%)\n', ...
    omp_range_rmse, 100*(rd_range_rmse - omp_range_rmse)/rd_range_rmse);
fprintf('组合估计: %.4f (vs 距离多普勒改进: %.1f%%)\n', ...
    combined_range_rmse, 100*(rd_range_rmse - combined_range_rmse)/rd_range_rmse);
fprintf('卡尔曼滤波: %.4f (vs 组合估计改进: %.1f%%)\n', ...
    kf_range_rmse, 100*(combined_range_rmse - kf_range_rmse)/combined_range_rmse);
fprintf('\n');

fprintf('方位角估计RMSE (度):\n');
fprintf('MUSIC: %.4f\n', music_azimuth_rmse);
fprintf('OMP稀疏重建: %.4f (vs MUSIC改进: %.1f%%)\n', ...
    omp_azimuth_rmse, 100*(music_azimuth_rmse - omp_azimuth_rmse)/music_azimuth_rmse);
fprintf('组合估计: %.4f (vs MUSIC改进: %.1f%%)\n', ...
    combined_azimuth_rmse, 100*(music_azimuth_rmse - combined_azimuth_rmse)/music_azimuth_rmse);
fprintf('卡尔曼滤波: %.4f (vs 组合估计改进: %.1f%%)\n', ...
    kf_azimuth_rmse, 100*(combined_azimuth_rmse - kf_azimuth_rmse)/combined_azimuth_rmse);
fprintf('\n');

fprintf('俯仰角估计RMSE (度):\n');
fprintf('MUSIC: %.4f\n', music_elevation_rmse);
fprintf('OMP稀疏重建: %.4f (vs MUSIC改进: %.1f%%)\n', ...
    omp_elevation_rmse, 100*(music_elevation_rmse - omp_elevation_rmse)/music_elevation_rmse);
fprintf('组合估计: %.4f (vs MUSIC改进: %.1f%%)\n', ...
    combined_elevation_rmse, 100*(music_elevation_rmse - combined_elevation_rmse)/music_elevation_rmse);
fprintf('卡尔曼滤波: %.4f (vs 组合估计改进: %.1f%%)\n', ...
    kf_elevation_rmse, 100*(combined_elevation_rmse - kf_elevation_rmse)/combined_elevation_rmse);
fprintf('\n');

fprintf('3D位置RMSE (米):\n');
fprintf('距离多普勒: %.4f\n', rd_3d_rmse);
fprintf('MUSIC: %.4f\n', music_3d_rmse);
fprintf('OMP稀疏重建: %.4f (vs 组合算法: %.1f%%)\n', ...
    omp_3d_rmse, 100*(combined_3d_rmse - omp_3d_rmse)/combined_3d_rmse);
fprintf('组合估计: %.4f\n', combined_3d_rmse);
fprintf('卡尔曼滤波: %.4f (vs 组合估计改进: %.1f%%)\n', ...
    kf_3d_rmse, 100*(combined_3d_rmse - kf_3d_rmse)/combined_3d_rmse);

% 可视化结果
try
    % 使用过渡帧数进行可视化
    visualize_results(true_positions, rd_positions, music_positions, combined_positions, omp_positions, estimated_positions, warmup_frames + smooth_frames);
    
    % 绘制卡尔曼滤波状态和协方差随时间变化
    if exist('state_history', 'var') && ~isempty(state_history)
        figure('Position', [100, 100, 1000, 600]);
        
        % 状态随时间变化
        subplot(2, 2, 1);
        plot(1:size(state_history, 1), state_history(:, 1), 'r-', 'LineWidth', 1.5);
        hold on;
        plot(1:size(state_history, 1), state_history(:, 4), 'b--', 'LineWidth', 1.5);
        title('距离和距离变化率');
        xlabel('帧');
        ylabel('数值');
        legend('距离 (m)', '距离变化率 (m/s)', 'Location', 'best');
        grid on;
        
        subplot(2, 2, 2);
        plot(1:size(state_history, 1), state_history(:, 2), 'r-', 'LineWidth', 1.5);
        hold on;
        plot(1:size(state_history, 1), state_history(:, 5), 'b--', 'LineWidth', 1.5);
        title('方位角和方位角变化率');
        xlabel('帧');
        ylabel('数值');
        legend('方位角 (度)', '方位角变化率 (度/s)', 'Location', 'best');
        grid on;
        
        subplot(2, 2, 3);
        plot(1:size(state_history, 1), state_history(:, 3), 'r-', 'LineWidth', 1.5);
        hold on;
        plot(1:size(state_history, 1), state_history(:, 6), 'b--', 'LineWidth', 1.5);
        title('俯仰角和俯仰角变化率');
        xlabel('帧');
        ylabel('数值');
        legend('俯仰角 (度)', '俯仰角变化率 (度/s)', 'Location', 'best');
        grid on;
        
        % 协方差矩阵对角线元素随时间变化
        if exist('covariance_history', 'var') && ~isempty(covariance_history)
            subplot(2, 2, 4);
            cov_diag = zeros(size(covariance_history, 1), 6);
            for i = 1:size(covariance_history, 1)
                cov_diag(i, :) = sqrt(diag(covariance_history{i}));
            end
            
            semilogy(1:size(cov_diag, 1), cov_diag(:, 1), 'r-', 'LineWidth', 1.5);
            hold on;
            semilogy(1:size(cov_diag, 1), cov_diag(:, 2), 'g--', 'LineWidth', 1.5);
            semilogy(1:size(cov_diag, 1), cov_diag(:, 3), 'b-.', 'LineWidth', 1.5);
            title('状态估计不确定性(标准差)');
            xlabel('帧');
            ylabel('标准差 (对数尺度)');
            legend('距离', '方位角', '俯仰角', 'Location', 'best');
            grid on;
        end
        
        sgtitle('卡尔曼滤波状态估计', 'FontSize', 14);
    end
catch e
    %warning('可视化结果时出错: %s', e.message);
    fprintf('尝试使用简化版可视化...\n');
    
    % 简单的可视化
    figure;
    
    % 绘制距离估计
    subplot(3, 1, 1);
    plot(1:size(true_positions, 1), true_positions(:, 1), 'k-', 'LineWidth', 2);
    hold on;
    plot(1:size(rd_positions, 1), rd_positions(:, 1), 'g:');
    plot(1:size(omp_positions, 1), omp_positions(:, 1), 'b--');
    plot(1:size(estimated_positions, 1), estimated_positions(:, 1), 'r-');
    title('距离估计');
    xlabel('帧');
    ylabel('距离 (m)');
    legend('真实值', '距离多普勒', 'OMP', '卡尔曼', 'Location', 'best');
    grid on;
    
    % 绘制方位角估计
    subplot(3, 1, 2);
    plot(1:size(true_positions, 1), true_positions(:, 2), 'k-', 'LineWidth', 2);
    hold on;
    plot(1:size(music_positions, 1), music_positions(:, 2), 'g:');
    plot(1:size(omp_positions, 1), omp_positions(:, 2), 'b--');
    plot(1:size(estimated_positions, 1), estimated_positions(:, 2), 'r-');
    title('方位角估计');
    xlabel('帧');
    ylabel('方位角 (度)');
    legend('真实值', 'MUSIC', 'OMP', '卡尔曼', 'Location', 'best');
    grid on;
    
    % 绘制俯仰角估计
    subplot(3, 1, 3);
    plot(1:size(true_positions, 1), true_positions(:, 3), 'k-', 'LineWidth', 2);
    hold on;
    plot(1:size(music_positions, 1), music_positions(:, 3), 'g:');
    plot(1:size(omp_positions, 1), omp_positions(:, 3), 'b--');
    plot(1:size(estimated_positions, 1), estimated_positions(:, 3), 'r-');
    title('俯仰角估计');
    xlabel('帧');
    ylabel('俯仰角 (度)');
    legend('真实值', 'MUSIC', 'OMP', '卡尔曼', 'Location', 'best');
    grid on;
end

function cart_coords = polar_to_cart(polar_coords)
% 将极坐标 [R, theta, phi] 转换为笛卡尔坐标 [x, y, z]
% R: 距离，theta: 方位角(度)，phi: 俯仰角(度)

n_points = size(polar_coords, 1);
cart_coords = zeros(n_points, 3);

for i = 1:n_points
    R = polar_coords(i, 1);
    theta = polar_coords(i, 2) * pi/180; % 转为弧度
    phi = polar_coords(i, 3) * pi/180;   % 转为弧度
    
    % 球坐标转笛卡尔坐标
    x = R * cos(phi) * cos(theta);
    y = R * cos(phi) * sin(theta);
    z = R * sin(phi);
    
    cart_coords(i, :) = [x, y, z];
end
end