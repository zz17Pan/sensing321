function atom = generate_atom_with_doppler(R, theta, phi, f_d, params, downsample_factor, residual)
% GENERATE_ATOM_WITH_DOPPLER 生成包含多普勒效应的字典原子
%   atom = GENERATE_ATOM_WITH_DOPPLER(R, theta, phi, f_d, params, downsample_factor, residual)
%   生成给定参数下的字典原子，包含多普勒频移
%   R: 距离(m)
%   theta: 方位角(弧度)
%   phi: 俯仰角(弧度)
%   f_d: 多普勒频移(Hz)
%   downsample_factor: 降采样因子
%   residual: 当前残差向量，用于维度匹配

% 使用单个接收子阵
n_rx_antennas = params.N_antennas_per_subarray;
n_chirps = params.N_chirps;
n_samples = params.N_samples;

% 降采样后的维度
n_chirps_ds = ceil(n_chirps / downsample_factor);
n_samples_ds = ceil(n_samples / downsample_factor);

% 计算原始分辨率下的时间点
t_chirp = (0:n_samples-1) / params.fs;

% 计算时延
tau = 2 * R / params.c;

% 方向矢量
direction = [cos(phi)*sin(theta), cos(phi)*cos(theta), sin(phi)];

% 初始化原子
atom = zeros(n_rx_antennas * n_chirps_ds * n_samples_ds, 1);

try
    % 获取接收子阵的天线位置
    rx_subarray_idx = params.sensing_rx_subarray;
    
    % 创建接收天线位置数组
    rx_array_pos = zeros(n_rx_antennas, 3);
    
    % 从接收子阵中提取天线位置
    ant_idx = 1;
    for nx = 1:4  % 假设4x4子阵
        for nz = 1:4
            % 计算相对于子阵中心的天线位置
            ant_offset = [(nx - 2.5) * params.d, 0, (nz - 2.5) * params.d];
            
            % 保存天线位置
            rx_array_pos(ant_idx, :) = ant_offset;
            ant_idx = ant_idx + 1;
        end
    end
    
    % 计算接收信号
    for i_chirp = 1:n_chirps_ds
        % 原始chirp索引
        orig_chirp = min((i_chirp-1)*downsample_factor + 1, n_chirps);
        
        % 计算chirp时间
        % 使用params.T（扫频时间）作为chirp时间
        if isfield(params, 'T_chirp')
            t_chirp_start = (orig_chirp-1) * params.T_chirp;
        else
            t_chirp_start = (orig_chirp-1) * params.T; % 使用扫频时间作为替代
        end
        
        for i_sample = 1:n_samples_ds
            % 原始sample索引
            orig_sample = min((i_sample-1)*downsample_factor + 1, n_samples);
            
            % 计算采样时间
            t = t_chirp_start + t_chirp(orig_sample);
            
            % 计算接收信号相位
            % 包含距离引起的时延和多普勒频移
            phase = -2*pi * (params.f0 * tau + 0.5 * params.k * tau^2 + f_d * t);
            
            % 计算接收信号
            for i_ant = 1:n_rx_antennas
                % 获取天线位置
                ant_pos = rx_array_pos(i_ant, :);
                
                % 计算方向引起的相位差
                phase_diff = 2*pi * dot(ant_pos, direction) / params.lambda;
                
                % 计算总相位
                total_phase = phase - phase_diff;
                
                % 计算接收信号
                idx = (i_chirp-1)*n_rx_antennas*n_samples_ds + (i_sample-1)*n_rx_antennas + i_ant;
                atom(idx) = exp(1j * total_phase);
            end
        end
    end
    
    % 归一化原子
    if norm(atom) > 0
        atom = atom / norm(atom);
    else
        fprintf('警告: 原子范数为零，无法归一化\n');
        % 创建一个随机单位向量作为替代
        atom = randn(size(atom)) + 1j*randn(size(atom));
        atom = atom / norm(atom);
    end
    
    % 确保原子维度与残差匹配
    if length(atom) ~= length(residual)
        fprintf('警告: 原子维度(%d)与残差维度(%d)不匹配，调整维度\n', length(atom), length(residual));
        if length(atom) > length(residual)
            % 截断原子
            atom = atom(1:length(residual));
        else
            % 扩展原子
            temp_atom = zeros(length(residual), 1);
            temp_atom(1:length(atom)) = atom;
            atom = temp_atom;
        end
        
        % 重新归一化
        if norm(atom) > 0
            atom = atom / norm(atom);
        end
    end
catch e
    % 错误处理
    fprintf('生成原子时出错: %s\n', e.message);
    % 创建一个与残差维度匹配的随机单位向量作为替代
    atom = randn(length(residual), 1) + 1j*randn(length(residual), 1);
    atom = atom / norm(atom);
end
end