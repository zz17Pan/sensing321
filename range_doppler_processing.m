function [range_est, velocity_est] = range_doppler_processing(rx_signal, params)
% RANGE_DOPPLER_PROCESSING 进行距离-多普勒处理和CFAR检测
%   [range_est, velocity_est] = RANGE_DOPPLER_PROCESSING(rx_signal, params)
%   对接收信号进行2D-FFT处理，生成距离-多普勒谱，并通过CFAR检测目标

% 获取信号维度
[n_rx_antennas, n_chirps, n_samples] = size(rx_signal);

% 检查输入信号是否有效
if isempty(rx_signal) || any(isnan(rx_signal(:))) || any(isinf(rx_signal(:)))
    fprintf('警告: 输入信号包含NaN或Inf，或为空。返回默认估计值\n');
    range_est = params.initial_R;
    velocity_est = 0;
    return;
end

% 初始化距离-多普勒谱
range_doppler_map = zeros(n_samples, n_chirps);

% 对子阵内所有天线的信号求平均 (非相干积累)
for ant_idx = 1:n_rx_antennas
    % 提取当前天线信号
    try
        % 重新排列为[样本,chirp]格式
        ant_signal = squeeze(rx_signal(ant_idx, :, :))';
        
        % 应用窗函数 (减少旁瓣)
        range_window = hamming(n_samples);
        doppler_window = hamming(n_chirps);
        
        windowed_signal = ant_signal .* range_window;
        
        % 距离FFT (快时间FFT)
        range_fft = fft(windowed_signal, [], 1);
        
        % 对每个距离bin应用多普勒窗口
        for r_bin = 1:n_samples
            range_fft(r_bin, :) = range_fft(r_bin, :) .* doppler_window';
        end
        
        % 多普勒FFT (慢时间FFT)
        range_doppler = fft(range_fft, [], 2);
        
        % 非相干积累
        range_doppler_map = range_doppler_map + abs(range_doppler).^2;
    catch e
        % 处理错误，继续执行
        fprintf('处理天线%d信号时出错: %s\n', ant_idx, e.message);
        continue;
    end
end

% 确保至少处理了一个天线信号
if all(range_doppler_map(:) == 0)
    fprintf('警告: 所有天线信号处理失败。返回默认估计值\n');
    range_est = params.initial_R;
    velocity_est = 0;
    return;
end

% 平均化
range_doppler_map = range_doppler_map / max(1, n_rx_antennas);

% 检查是否包含NaN或Inf
if any(isnan(range_doppler_map(:))) || any(isinf(range_doppler_map(:)))
    fprintf('警告: 距离-多普勒谱包含NaN或Inf，将替换为0\n');
    range_doppler_map(isnan(range_doppler_map) | isinf(range_doppler_map)) = 0;
end

% 调整频谱，使零多普勒位于中心
range_doppler_map = fftshift(range_doppler_map, 2);

% CFAR检测
try
    [peaks, peak_indices] = cfar_detector(range_doppler_map, params);
catch e
    fprintf('CFAR检测出错: %s\n', e.message);
    % 在CFAR失败时，直接寻找最强点
    peaks = [];
    peak_indices = [];
end

% 如果未检测到峰值或CFAR失败，选择最强点
if isempty(peaks) || isempty(peak_indices)
    fprintf('使用最强点作为目标\n');
    [~, max_idx] = max(range_doppler_map(:));
    [r_bin, d_bin] = ind2sub(size(range_doppler_map), max_idx);
    peak_indices = [r_bin, d_bin];
    peaks = range_doppler_map(r_bin, d_bin);
end

% 取出最强峰值
[~, max_peak_idx] = max(peaks);
r_bin = peak_indices(max_peak_idx, 1);
d_bin = peak_indices(max_peak_idx, 2);

% 计算估计的距离
freq_res = params.fs / n_samples;
range_res = params.c / (2 * params.B);

% 确保r_bin在有效范围内
if r_bin > n_samples/2
    r_bin = r_bin - n_samples; % 负频率校正
end
range_est = max(0.1, abs(r_bin) * range_res); % 确保距离为正值且不为零

% 计算估计的速度
doppler_res = 1 / (n_chirps * params.T);
velocity_res = params.lambda / 2 * doppler_res;
d_bin_centered = d_bin - n_chirps/2 - 1; % 中心化
velocity_est = d_bin_centered * velocity_res;

% 打印估计结果
fprintf('距离-多普勒处理: 距离估计=%.2fm, 速度估计=%.2fm/s\n', range_est, velocity_est);

end

function [peaks, peak_indices] = cfar_detector(range_doppler_map, params)
% 使用CA-CFAR (Cell-Averaging CFAR) 检测目标
% 参数解包
Pfa = params.Pfa;
guard_cells = params.guard_cells;
training_cells = params.training_cells;

% 获取数据尺寸
[n_range_bins, n_doppler_bins] = size(range_doppler_map);

% 初始化输出
peaks = [];
peak_indices = [];

% 计算CFAR窗口大小
range_guard = guard_cells(1);
doppler_guard = guard_cells(2);
range_train = training_cells(1);
doppler_train = training_cells(2);

% 检查CFAR窗口是否合理
if 2*(range_guard+range_train) >= n_range_bins || 2*(doppler_guard+doppler_train) >= n_doppler_bins
    fprintf('警告: CFAR窗口尺寸过大，无法进行CFAR检测\n');
    % 使用最大值替代CFAR
    [~, max_idx] = max(range_doppler_map(:));
    [r_idx, d_idx] = ind2sub(size(range_doppler_map), max_idx);
    peaks = range_doppler_map(r_idx, d_idx);
    peak_indices = [r_idx, d_idx];
    return;
end

% 设置从训练单元数量计算的CFAR缩放因子
% 根据恒虚警率的要求，我们计算缩放因子
num_training_cells = (2*range_train + 2*doppler_train + 4*range_train*doppler_train);
alpha = max(1.0, num_training_cells * (Pfa^(-1/num_training_cells) - 1)); % 确保alpha不小于1.0

% 对每个单元执行CFAR处理
for r_idx = 1+range_guard+range_train : n_range_bins-range_guard-range_train
    for d_idx = 1+doppler_guard+doppler_train : n_doppler_bins-doppler_guard-doppler_train
        % 当前单元值
        cell_under_test = range_doppler_map(r_idx, d_idx);
        
        % 检查单元格是否有效
        if isnan(cell_under_test) || isinf(cell_under_test)
            continue;
        end
        
        % 提取训练单元区域 (排除保护单元)
        % 确保所有提取区域都在有效范围内
        
        % 左侧训练区域
        left_train_r_start = max(1, r_idx-range_guard-range_train);
        left_train_r_end = max(1, r_idx-range_guard-1);
        left_train_d_start = max(1, d_idx-doppler_guard-doppler_train);
        left_train_d_end = min(n_doppler_bins, d_idx+doppler_guard+doppler_train);
        
        left_train = range_doppler_map(left_train_r_start:left_train_r_end, ...
                                       left_train_d_start:left_train_d_end);
        
        % 右侧训练区域
        right_train_r_start = min(n_range_bins, r_idx+range_guard+1);
        right_train_r_end = min(n_range_bins, r_idx+range_guard+range_train);
        right_train_d_start = max(1, d_idx-doppler_guard-doppler_train);
        right_train_d_end = min(n_doppler_bins, d_idx+doppler_guard+doppler_train);
        
        right_train = range_doppler_map(right_train_r_start:right_train_r_end, ...
                                        right_train_d_start:right_train_d_end);
        
        % 上方训练区域
        top_train_r_start = max(1, r_idx-range_guard);
        top_train_r_end = min(n_range_bins, r_idx+range_guard);
        top_train_d_start = min(n_doppler_bins, d_idx+doppler_guard+1);
        top_train_d_end = min(n_doppler_bins, d_idx+doppler_guard+doppler_train);
        
        top_train = range_doppler_map(top_train_r_start:top_train_r_end, ...
                                      top_train_d_start:top_train_d_end);
        
        % 下方训练区域
        bottom_train_r_start = max(1, r_idx-range_guard);
        bottom_train_r_end = min(n_range_bins, r_idx+range_guard);
        bottom_train_d_start = max(1, d_idx-doppler_guard-doppler_train);
        bottom_train_d_end = max(1, d_idx-doppler_guard-1);
        
        bottom_train = range_doppler_map(bottom_train_r_start:bottom_train_r_end, ...
                                         bottom_train_d_start:bottom_train_d_end);
        
        % 合并所有训练区域
        all_train_cells = [left_train(:); right_train(:); top_train(:); bottom_train(:)];
        
        % 移除NaN和Inf
        all_train_cells = all_train_cells(~isnan(all_train_cells) & ~isinf(all_train_cells));
        
        % 确保有足够的训练单元
        if length(all_train_cells) < 0.5 * num_training_cells
            continue;
        end
        
        % 计算所有训练单元的平均值
        noise_level = mean(all_train_cells) + eps; % 添加eps避免除零
        
        % 设置自适应阈值
        threshold = alpha * noise_level;
        
        % 检测是否超过阈值
        if cell_under_test > threshold
            % 局部最大值检测 (确保是峰值)
            local_r_start = max(1, r_idx-1);
            local_r_end = min(n_range_bins, r_idx+1);
            local_d_start = max(1, d_idx-1);
            local_d_end = min(n_doppler_bins, d_idx+1);
            
            local_region = range_doppler_map(local_r_start:local_r_end, ...
                                             local_d_start:local_d_end);
            
            if cell_under_test == max(local_region(:))
                peaks = [peaks; cell_under_test];
                peak_indices = [peak_indices; r_idx, d_idx];
            end
        end
    end
end

% 确保至少返回一个峰值
if isempty(peaks)
    fprintf('CFAR未检测到目标，使用全局最大值\n');
    [~, max_idx] = max(range_doppler_map(:));
    [r_idx, d_idx] = ind2sub(size(range_doppler_map), max_idx);
    peaks = range_doppler_map(r_idx, d_idx);
    peak_indices = [r_idx, d_idx];
end

end 