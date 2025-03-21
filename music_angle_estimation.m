function [theta_est, phi_est] = music_angle_estimation(rx_signal, params)
% MUSIC_ANGLE_ESTIMATION 使用MUSIC算法估计目标的方位角和俯仰角
%   [theta_est, phi_est] = MUSIC_ANGLE_ESTIMATION(rx_signal, params)
%   对接收信号进行MUSIC算法处理，估计目标的方位角和俯仰角

% 获取信号维度
[n_rx_antennas, n_chirps, n_samples] = size(rx_signal);

% 对每个距离bin执行距离FFT
range_fft = zeros(n_rx_antennas, n_chirps, n_samples);
for ant_idx = 1:n_rx_antennas
    ant_signal = squeeze(rx_signal(ant_idx, :, :));
    range_fft(ant_idx, :, :) = fft(ant_signal, [], 2);
end

% 计算距离谱
range_spectrum = mean(abs(range_fft).^2, 2);
range_spectrum = squeeze(mean(range_spectrum, 1));

% 找到距离谱中的最强峰值对应的距离bin
[~, max_range_bin] = max(range_spectrum);

% 应用平滑处理，减少噪声影响
range_spectrum_smooth = conv(range_spectrum, ones(5,1)/5, 'same');
[~, max_range_bin_smooth] = max(range_spectrum_smooth);

% 使用平滑后的峰值位置
max_range_bin = max_range_bin_smooth;

% 从该距离bin提取信号
signals_at_peak = squeeze(range_fft(:, :, max_range_bin));

% 估计协方差矩阵
R = signals_at_peak * signals_at_peak' / n_chirps;

% 对协方差矩阵进行特征分解
[V, D] = eig(R);
eigenvalues = diag(D);

% 按特征值降序排序
[eigenvalues, idx] = sort(eigenvalues, 'descend');
V = V(:, idx);

% 动态估计信号源数量
% 使用MDL(最小描述长度)准则估计信号源数量
[n_sources_mdl] = estimate_source_number(eigenvalues, n_chirps);

% 设置最大可能的信号源数量为天线数量的一半
max_sources = floor(n_rx_antennas/2);
n_sources = min(max_sources, max(1, n_sources_mdl));

fprintf('MUSIC: 估计信号源数量 = %d\n', n_sources);

% 区分信号子空间和噪声子空间
noise_subspace = V(:, n_sources+1:end);

% 创建角度搜索网格
theta_range = params.theta_range;
phi_range = params.phi_range;
theta_grid = theta_range(1):params.angle_grid_step:theta_range(2);
phi_grid = phi_range(1):params.angle_grid_step:phi_range(2);

% 初始化MUSIC谱
n_theta = length(theta_grid);
n_phi = length(phi_grid);
music_spectrum = zeros(n_theta, n_phi);

% 计算MUSIC谱
for i = 1:n_theta
    theta = theta_grid(i) * pi/180;  % 转换为弧度
    
    for j = 1:n_phi
        phi = phi_grid(j) * pi/180;  % 转换为弧度
        
        % 计算阵列响应向量 (仅用于单个子阵)
        a = array_response_vector(theta, phi, params);
        
        % 计算MUSIC谱
        proj = noise_subspace' * a;
        music_spectrum(i, j) = 1 / (proj' * proj + eps);  % 添加eps避免除零
    end
end

% 找到MUSIC谱中的峰值
[~, max_idx] = max(music_spectrum(:));
[theta_idx, phi_idx] = ind2sub([n_theta, n_phi], max_idx);

% 转换为角度估计
theta_est = theta_grid(theta_idx);
phi_est = phi_grid(phi_idx);

end

function a = array_response_vector(theta, phi, params)
% 计算阵列响应向量 (仅用于单个子阵)
% theta: 方位角(弧度)
% phi: 俯仰角(弧度)

% 子阵内天线数量
Nx = 4;  % 子阵x方向天线数
Nz = 4;  % 子阵z方向天线数
n_antennas = Nx * Nz;

% 初始化阵列响应向量
a = zeros(n_antennas, 1);

% 计算方向余弦
u = cos(phi) * sin(theta);
v = sin(phi);
w = cos(phi) * cos(theta);

% 生成子阵内每个天线的响应
ant_idx = 1;
for nz = 1:Nz
    for nx = 1:Nx
        % 天线相对于子阵中心的位置
        x = (nx - 2.5) * params.d;
        z = (nz - 2.5) * params.d;
        
        % 计算相位
        phase = 2 * pi / params.lambda * (x * u + z * w);
        
        % 计算阵列响应
        a(ant_idx) = exp(1j * phase);
        
        ant_idx = ant_idx + 1;
    end
end

% 归一化
a = a / norm(a);

end

function [n_sources] = estimate_source_number(eigenvalues, n_samples)
% 使用MDL(最小描述长度)准则估计信号源数量
% eigenvalues: 协方差矩阵的特征值
% n_samples: 样本数量

L = length(eigenvalues);
mdl_values = zeros(L-1, 1);

for k = 0:(L-2)
    % 分离噪声特征值
    lambda_noise = eigenvalues(k+1:end);
    
    % 计算MDL
    if k == 0
        G_k = 0;
    else
        G_k = -n_samples * sum(log(lambda_noise / (sum(lambda_noise) / length(lambda_noise))));
    end
    
    % 自由参数数量
    P_k = k * (2*L - k);
    
    % MDL值
    mdl_values(k+1) = G_k + 0.5 * P_k * log(n_samples);
end

% 找到最小MDL值对应的信号源数量
[~, min_idx] = min(mdl_values);
n_sources = min_idx;

end