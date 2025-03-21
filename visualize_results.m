function visualize_results(true_positions, rd_positions, music_positions, combined_positions, omp_positions, estimated_positions, transition_frames)
% VISUALIZE_RESULTS 可视化真实轨迹和各种估计轨迹
%   VISUALIZE_RESULTS(true_positions, rd_positions, music_positions, combined_positions, omp_positions, estimated_positions, transition_frames) 
%   绘制真实位置和各种算法估计位置随时间的变化，并标记预热和平滑过渡阶段

% 检查输入参数
if nargin < 7
    transition_frames = 10; % 默认过渡帧数
    fprintf('未提供过渡帧数，默认为%d\n', transition_frames);
end

% 检查输入数据的有效性
if isempty(true_positions) || isempty(estimated_positions)
    error('输入轨迹数据为空');
end

% 获取轨迹长度
n_frames = size(true_positions, 1);
time_axis = 1:n_frames;

% 检查所有输入数据维度是否一致
if size(rd_positions, 1) ~= n_frames || size(music_positions, 1) ~= n_frames || ...
   size(combined_positions, 1) ~= n_frames || size(omp_positions, 1) ~= n_frames || ...
   size(estimated_positions, 1) ~= n_frames
    warning('输入数据维度不一致！调整为相同维度');
    % 调整为最小长度
    min_len = min([size(true_positions, 1), size(rd_positions, 1), ...
                 size(music_positions, 1), size(combined_positions, 1), ...
                 size(omp_positions, 1), size(estimated_positions, 1)]);
    true_positions = true_positions(1:min_len, :);
    rd_positions = rd_positions(1:min_len, :);
    music_positions = music_positions(1:min_len, :);
    combined_positions = combined_positions(1:min_len, :);
    omp_positions = omp_positions(1:min_len, :);
    estimated_positions = estimated_positions(1:min_len, :);
    n_frames = min_len;
    time_axis = 1:n_frames;
end

% 确保transition_frames不超过总帧数的一半
transition_frames = min(transition_frames, floor(n_frames/2));

% 创建图形
try
    figure('Position', [100, 100, 1200, 800]);
    
    % 1. 距离随时间变化
    subplot(3, 2, 1);
    plot(time_axis, true_positions(:, 1), 'k-', 'LineWidth', 2);
    hold on;
    plot(time_axis, rd_positions(:, 1), 'g:', 'LineWidth', 1.5);
    plot(time_axis, combined_positions(:, 1), 'm-.', 'LineWidth', 1.5);
    plot(time_axis, omp_positions(:, 1), 'b--', 'LineWidth', 1.5);
    plot(time_axis, estimated_positions(:, 1), 'r-', 'LineWidth', 1.5);
    
    % 标记预热和过渡阶段
    xline(transition_frames + 0.5, '--', '过渡结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    fill([0 transition_frames+0.5 transition_frames+0.5 0], [min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    grid on;
    xlabel('帧');
    ylabel('距离 (m)');
    title('距离随时间变化');
    legend('真实值', '距离多普勒', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');
    
    % 2. 距离误差随时间变化
    subplot(3, 2, 2);
    rd_range_error = abs(rd_positions(:, 1) - true_positions(:, 1));
    combined_range_error = abs(combined_positions(:, 1) - true_positions(:, 1));
    omp_range_error = abs(omp_positions(:, 1) - true_positions(:, 1));
    kf_range_error = abs(estimated_positions(:, 1) - true_positions(:, 1));
    
    plot(time_axis, rd_range_error, 'g:', 'LineWidth', 1.5);
    hold on;
    plot(time_axis, combined_range_error, 'm-.', 'LineWidth', 1.5);
    plot(time_axis, omp_range_error, 'b--', 'LineWidth', 1.5);
    plot(time_axis, kf_range_error, 'r-', 'LineWidth', 1.5);
    
    % 标记预热和过渡阶段
    xline(transition_frames + 0.5, '--', '过渡结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    fill([0 transition_frames+0.5 transition_frames+0.5 0], [min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    grid on;
    xlabel('帧');
    ylabel('距离误差 (m)');
    title('距离估计误差(绝对值)');
    legend('距离多普勒', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');
    
    % 3. 方位角随时间变化
    subplot(3, 2, 3);
    plot(time_axis, true_positions(:, 2), 'k-', 'LineWidth', 2);
    hold on;
    plot(time_axis, music_positions(:, 2), 'g:', 'LineWidth', 1.5);
    plot(time_axis, combined_positions(:, 2), 'm-.', 'LineWidth', 1.5);
    plot(time_axis, omp_positions(:, 2), 'b--', 'LineWidth', 1.5);
    plot(time_axis, estimated_positions(:, 2), 'r-', 'LineWidth', 1.5);
    
    % 标记预热和过渡阶段
    xline(transition_frames + 0.5, '--', '过渡结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    fill([0 transition_frames+0.5 transition_frames+0.5 0], [min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    grid on;
    xlabel('帧');
    ylabel('方位角 (度)');
    title('方位角随时间变化');
    legend('真实值', 'MUSIC', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');
    
    % 4. 方位角误差随时间变化
    subplot(3, 2, 4);
    music_azimuth_error = abs(music_positions(:, 2) - true_positions(:, 2));
    combined_azimuth_error = abs(combined_positions(:, 2) - true_positions(:, 2));
    omp_azimuth_error = abs(omp_positions(:, 2) - true_positions(:, 2));
    kf_azimuth_error = abs(estimated_positions(:, 2) - true_positions(:, 2));
    
    plot(time_axis, music_azimuth_error, 'g:', 'LineWidth', 1.5);
    hold on;
    plot(time_axis, combined_azimuth_error, 'm-.', 'LineWidth', 1.5);
    plot(time_axis, omp_azimuth_error, 'b--', 'LineWidth', 1.5);
    plot(time_axis, kf_azimuth_error, 'r-', 'LineWidth', 1.5);
    
    % 标记预热和过渡阶段
    xline(transition_frames + 0.5, '--', '过渡结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    fill([0 transition_frames+0.5 transition_frames+0.5 0], [min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    grid on;
    xlabel('帧');
    ylabel('方位角误差 (度)');
    title('方位角估计误差(绝对值)');
    legend('MUSIC', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');
    
    % 5. 俯仰角随时间变化
    subplot(3, 2, 5);
    plot(time_axis, true_positions(:, 3), 'k-', 'LineWidth', 2);
    hold on;
    plot(time_axis, music_positions(:, 3), 'g:', 'LineWidth', 1.5);
    plot(time_axis, combined_positions(:, 3), 'm-.', 'LineWidth', 1.5);
    plot(time_axis, omp_positions(:, 3), 'b--', 'LineWidth', 1.5);
    plot(time_axis, estimated_positions(:, 3), 'r-', 'LineWidth', 1.5);
    
    % 标记预热和过渡阶段
    xline(transition_frames + 0.5, '--', '过渡结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    fill([0 transition_frames+0.5 transition_frames+0.5 0], [min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    grid on;
    xlabel('帧');
    ylabel('俯仰角 (度)');
    title('俯仰角随时间变化');
    legend('真实值', 'MUSIC', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');
    
    % 6. 俯仰角误差随时间变化
    subplot(3, 2, 6);
    music_elevation_error = abs(music_positions(:, 3) - true_positions(:, 3));
    combined_elevation_error = abs(combined_positions(:, 3) - true_positions(:, 3));
    omp_elevation_error = abs(omp_positions(:, 3) - true_positions(:, 3));
    kf_elevation_error = abs(estimated_positions(:, 3) - true_positions(:, 3));
    
    plot(time_axis, music_elevation_error, 'g:', 'LineWidth', 1.5);
    hold on;
    plot(time_axis, combined_elevation_error, 'm-.', 'LineWidth', 1.5);
    plot(time_axis, omp_elevation_error, 'b--', 'LineWidth', 1.5);
    plot(time_axis, kf_elevation_error, 'r-', 'LineWidth', 1.5);
    
    % 标记预热和过渡阶段
    xline(transition_frames + 0.5, '--', '过渡结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    fill([0 transition_frames+0.5 transition_frames+0.5 0], [min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    grid on;
    xlabel('帧');
    ylabel('俯仰角误差 (度)');
    title('俯仰角估计误差(绝对值)');
    legend('MUSIC', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');
    
    % 调整整体布局
    sgtitle('感知辅助太赫兹波束对准：参数估计性能', 'FontSize', 16);
    
    % 创建3D轨迹图
    try
        figure('Position', [100, 100, 900, 700]);
        
        % 将极坐标转换为笛卡尔坐标
        true_cart = polar_to_cartesian(true_positions);
        kf_cart = polar_to_cartesian(estimated_positions);
        omp_cart = polar_to_cartesian(omp_positions);
        combined_cart = polar_to_cartesian(combined_positions);
        
        % 绘制3D轨迹
        plot3(true_cart(:,1), true_cart(:,2), true_cart(:,3), 'k-', 'LineWidth', 2);
        hold on;
        plot3(combined_cart(:,1), combined_cart(:,2), combined_cart(:,3), 'm-.', 'LineWidth', 1.5);
        plot3(omp_cart(:,1), omp_cart(:,2), omp_cart(:,3), 'b--', 'LineWidth', 1.5);
        plot3(kf_cart(:,1), kf_cart(:,2), kf_cart(:,3), 'r-', 'LineWidth', 1.5);
        
        % 标记起始点和结束点
        plot3(true_cart(1,1), true_cart(1,2), true_cart(1,3), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
        plot3(true_cart(end,1), true_cart(end,2), true_cart(end,3), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        text(true_cart(1,1), true_cart(1,2), true_cart(1,3), '  起点', 'FontSize', 12);
        text(true_cart(end,1), true_cart(end,2), true_cart(end,3), '  终点', 'FontSize', 12);
        
        % 标记过渡阶段结束点
        if transition_frames < n_frames
            plot3(true_cart(transition_frames,1), true_cart(transition_frames,2), true_cart(transition_frames,3), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
            text(true_cart(transition_frames,1), true_cart(transition_frames,2), true_cart(transition_frames,3), '  过渡结束', 'FontSize', 12);
        end
        
        % 设置视角和网格
        grid on;
        xlabel('X (m)');
        ylabel('Y (m)');
        zlabel('Z (m)');
        title('3D轨迹可视化', 'FontSize', 14);
        legend('真实轨迹', '组合估计', 'OMP估计', '卡尔曼滤波', 'Location', 'best');
        view(30, 30);
        axis equal;
        
        % 计算RMSE - 仅使用过渡阶段后的数据
        eval_frames = (transition_frames+1):n_frames;
        
        if ~isempty(eval_frames) && length(eval_frames) > 1
            % 计算欧氏距离误差
            combined_dist_error = sqrt(sum((true_cart(eval_frames,:) - combined_cart(eval_frames,:)).^2, 2));
            omp_dist_error = sqrt(sum((true_cart(eval_frames,:) - omp_cart(eval_frames,:)).^2, 2));
            kf_dist_error = sqrt(sum((true_cart(eval_frames,:) - kf_cart(eval_frames,:)).^2, 2));
            
            % 计算RMSE
            combined_rmse = sqrt(mean(combined_dist_error.^2));
            omp_rmse = sqrt(mean(omp_dist_error.^2));
            kf_rmse = sqrt(mean(kf_dist_error.^2));
            
            % 显示RMSE
            text_x = min(get(gca, 'xlim')) + 0.1 * (max(get(gca, 'xlim')) - min(get(gca, 'xlim')));
            text_y = min(get(gca, 'ylim')) + 0.1 * (max(get(gca, 'ylim')) - min(get(gca, 'ylim')));
            text_z = max(get(gca, 'zlim'));
            
            text_str = sprintf('3D RMSE (过渡后): 组合=%.2fm, OMP=%.2fm, 卡尔曼=%.2fm', ...
                combined_rmse, omp_rmse, kf_rmse);
            text(text_x, text_y, text_z, text_str, 'FontSize', 12, 'Color', 'blue');
        end
        
    catch e
        warning('绘制3D轨迹图出错: %s', e.message);
    end
    
    % 绘制误差随时间变化的RMSE曲线 (使用滑动窗口)
    try
        figure('Position', [100, 100, 900, 400]);
        
        % 设置滑动窗口大小
        window_size = min(10, floor(n_frames/5));
        window_size = max(window_size, 2); % 确保窗口至少包含2个点
        
        % 初始化滑动RMSE数组
        sliding_frames = window_size:n_frames;
        rd_sliding_rmse = zeros(length(sliding_frames), 3);
        music_sliding_rmse = zeros(length(sliding_frames), 3);
        combined_sliding_rmse = zeros(length(sliding_frames), 3);
        omp_sliding_rmse = zeros(length(sliding_frames), 3);
        kf_sliding_rmse = zeros(length(sliding_frames), 3);
        
        % 计算滑动窗口RMSE
        for i = 1:length(sliding_frames)
            frame_idx = sliding_frames(i);
            window_indices = (frame_idx-window_size+1):frame_idx;
            
            % 距离RMSE
            rd_sliding_rmse(i, 1) = sqrt(mean((true_positions(window_indices, 1) - rd_positions(window_indices, 1)).^2));
            combined_sliding_rmse(i, 1) = sqrt(mean((true_positions(window_indices, 1) - combined_positions(window_indices, 1)).^2));
            omp_sliding_rmse(i, 1) = sqrt(mean((true_positions(window_indices, 1) - omp_positions(window_indices, 1)).^2));
            kf_sliding_rmse(i, 1) = sqrt(mean((true_positions(window_indices, 1) - estimated_positions(window_indices, 1)).^2));
            
            % 方位角RMSE
            music_sliding_rmse(i, 2) = sqrt(mean((true_positions(window_indices, 2) - music_positions(window_indices, 2)).^2));
            combined_sliding_rmse(i, 2) = sqrt(mean((true_positions(window_indices, 2) - combined_positions(window_indices, 2)).^2));
            omp_sliding_rmse(i, 2) = sqrt(mean((true_positions(window_indices, 2) - omp_positions(window_indices, 2)).^2));
            kf_sliding_rmse(i, 2) = sqrt(mean((true_positions(window_indices, 2) - estimated_positions(window_indices, 2)).^2));
            
            % 俯仰角RMSE
            music_sliding_rmse(i, 3) = sqrt(mean((true_positions(window_indices, 3) - music_positions(window_indices, 3)).^2));
            combined_sliding_rmse(i, 3) = sqrt(mean((true_positions(window_indices, 3) - combined_positions(window_indices, 3)).^2));
            omp_sliding_rmse(i, 3) = sqrt(mean((true_positions(window_indices, 3) - omp_positions(window_indices, 3)).^2));
            kf_sliding_rmse(i, 3) = sqrt(mean((true_positions(window_indices, 3) - estimated_positions(window_indices, 3)).^2));
        end
        
        % 绘制滑动RMSE曲线
        subplot(1, 3, 1);
        plot(sliding_frames, rd_sliding_rmse(:, 1), 'g:', 'LineWidth', 1.5);
        hold on;
        plot(sliding_frames, combined_sliding_rmse(:, 1), 'm-.', 'LineWidth', 1.5);
        plot(sliding_frames, omp_sliding_rmse(:, 1), 'b--', 'LineWidth', 1.5);
        plot(sliding_frames, kf_sliding_rmse(:, 1), 'r-', 'LineWidth', 1.5);
        
        % 标记过渡阶段结束
        if transition_frames < n_frames
            xline(transition_frames + 0.5, '--', '过渡结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
        end
        
        grid on;
        xlabel('帧');
        ylabel('滑动RMSE (m)');
        title('距离估计滑动RMSE');
        legend('距离多普勒', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');
        
        subplot(1, 3, 2);
        plot(sliding_frames, music_sliding_rmse(:, 2), 'g:', 'LineWidth', 1.5);
        hold on;
        plot(sliding_frames, combined_sliding_rmse(:, 2), 'm-.', 'LineWidth', 1.5);
        plot(sliding_frames, omp_sliding_rmse(:, 2), 'b--', 'LineWidth', 1.5);
        plot(sliding_frames, kf_sliding_rmse(:, 2), 'r-', 'LineWidth', 1.5);
        
        % 标记过渡阶段结束
        if transition_frames < n_frames
            xline(transition_frames + 0.5, '--', '过渡结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
        end
        
        grid on;
        xlabel('帧');
        ylabel('滑动RMSE (度)');
        title('方位角估计滑动RMSE');
        legend('MUSIC', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');
        
        subplot(1, 3, 3);
        plot(sliding_frames, music_sliding_rmse(:, 3), 'g:', 'LineWidth', 1.5);
        hold on;
        plot(sliding_frames, combined_sliding_rmse(:, 3), 'm-.', 'LineWidth', 1.5);
        plot(sliding_frames, omp_sliding_rmse(:, 3), 'b--', 'LineWidth', 1.5);
        plot(sliding_frames, kf_sliding_rmse(:, 3), 'r-', 'LineWidth', 1.5);
        
        % 标记过渡阶段结束
        if transition_frames < n_frames
            xline(transition_frames + 0.5, '--', '过渡结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
        end
        
        grid on;
        xlabel('帧');
        ylabel('滑动RMSE (度)');
        title('俯仰角估计滑动RMSE');
        legend('MUSIC', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');
        
        sgtitle(sprintf('参数估计滑动RMSE (窗口大小=%d帧)', window_size), 'FontSize', 14);
    catch e
        warning('绘制RMSE曲线出错: %s', e.message);
    end
    
catch e
    warning('可视化过程中出错: %s', e.message);
    % 创建简单的替代图形
    figure;
    
    % 绘制距离跟踪结果
    subplot(3, 1, 1);
    plot(time_axis, true_positions(:, 1), 'k-', 'LineWidth', 2);
    hold on;
    plot(time_axis, combined_positions(:, 1), 'g--');
    plot(time_axis, omp_positions(:, 1), 'b-.');
    plot(time_axis, estimated_positions(:, 1), 'r-');
    legend('真实', '组合', 'OMP', '卡尔曼');
    title('距离跟踪');
    xlabel('帧');
    ylabel('距离(m)');
    grid on;
    
    % 绘制方位角跟踪结果
    subplot(3, 1, 2);
    plot(time_axis, true_positions(:, 2), 'k-', 'LineWidth', 2);
    hold on;
    plot(time_axis, combined_positions(:, 2), 'g--');
    plot(time_axis, omp_positions(:, 2), 'b-.');
    plot(time_axis, estimated_positions(:, 2), 'r-');
    title('方位角跟踪');
    xlabel('帧');
    ylabel('方位角(度)');
    grid on;
    
    % 绘制俯仰角跟踪结果
    subplot(3, 1, 3);
    plot(time_axis, true_positions(:, 3), 'k-', 'LineWidth', 2);
    hold on;
    plot(time_axis, combined_positions(:, 3), 'g--');
    plot(time_axis, omp_positions(:, 3), 'b-.');
    plot(time_axis, estimated_positions(:, 3), 'r-');
    title('俯仰角跟踪');
    xlabel('帧');
    ylabel('俯仰角(度)');
    grid on;
    
    sgtitle('参数估计和跟踪性能（简化版）');
end
end

function cart_coords = polar_to_cartesian(polar_coords)
% 将极坐标转换为笛卡尔坐标
% polar_coords: [R, theta, phi] - R:距离，theta:方位角(度)，phi:俯仰角(度)

n_points = size(polar_coords, 1);
cart_coords = zeros(n_points, 3);

for i = 1:n_points
    R = polar_coords(i, 1);
    theta = polar_coords(i, 2) * pi/180; % 转换为弧度
    phi = polar_coords(i, 3) * pi/180;   % 转换为弧度
    
    % 转换公式: 
    % x = R * cos(phi) * cos(theta)
    % y = R * cos(phi) * sin(theta)
    % z = R * sin(phi)
    cart_coords(i, 1) = R * cos(phi) * cos(theta);
    cart_coords(i, 2) = R * cos(phi) * sin(theta);
    cart_coords(i, 3) = R * sin(phi);
end
end