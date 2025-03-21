function [R_est, theta_est, phi_est] = omp_sparse_reconstruction(rx_signal, prior_params, params)
% OMP_SPARSE_RECONSTRUCTION 使用先验信息的OMP算法进行参数估计
%   [R_est, theta_est, phi_est] = OMP_SPARSE_RECONSTRUCTION(rx_signal, prior_params, params)
%   利用距离-多普勒和MUSIC角度估计结果作为先验，通过OMP算法进行稀疏重建

% 获取先验信息
R_prior = prior_params.R_est;
theta_prior = prior_params.theta_est;
phi_prior = prior_params.phi_est;
sigma_R = prior_params.sigma_R;
sigma_theta = prior_params.sigma_theta;
sigma_phi = prior_params.sigma_phi;

% 定义多普勒速度假设选项 - 减少选项数量提高效率
v_r_options = [0.0, 5.0, 10.0, 15.0]; % 减少速度选项，节省计算资源
f_d_options = 2 * v_r_options * params.fc / params.c;

% ===== 优化搜索范围和分辨率，提高精度和效率 =====
% 扩大搜索范围 (以先验为中心，覆盖3个标准差)
R_range = [max(0.1, R_prior - 3*sigma_R), R_prior + 3*sigma_R];
theta_range = [theta_prior - 3*sigma_theta, theta_prior + 3*sigma_theta];
phi_range = [phi_prior - 3*sigma_phi, phi_prior + 3*sigma_phi];

% 打印搜索范围
fprintf('OMP搜索范围: 距离=[%.2f, %.2f]m, 方位角=[%.2f, %.2f]°, 俯仰角=[%.2f, %.2f]°\n', ...
    R_range(1), R_range(2), theta_range(1), theta_range(2), phi_range(1), phi_range(2));

% 确保角度范围在合理区间内
theta_range(1) = max(params.theta_range(1), theta_range(1));
theta_range(2) = min(params.theta_range(2), theta_range(2));
phi_range(1) = max(params.phi_range(1), phi_range(1));
phi_range(2) = min(params.phi_range(2), phi_range(2));

% 优化网格步长，平衡精度和计算量
% 根据搜索范围动态调整网格密度
R_range_size = R_range(2) - R_range(1);
R_grid_step = max(params.R_grid_step * 0.8, R_range_size / 10);

theta_range_size = theta_range(2) - theta_range(1);
theta_grid_step = max(params.angle_grid_step * 0.8, theta_range_size / 10);

phi_range_size = phi_range(2) - phi_range(1);
phi_grid_step = max(params.angle_grid_step * 0.8, phi_range_size / 10);

% 定义参数网格 - 针对先验值附近使用更密集的网格
R_grid = create_nonuniform_grid(R_range(1), R_range(2), R_prior, R_grid_step);
theta_grid = create_nonuniform_grid(theta_range(1), theta_range(2), theta_prior, theta_grid_step);
phi_grid = create_nonuniform_grid(phi_range(1), phi_range(2), phi_prior, phi_grid_step);

% 获取信号维度
[n_rx_antennas, n_chirps, n_samples] = size(rx_signal);

% 降低信号采样维度以提高效率，但不要降低得太多
downsample_factor = 2;  % 降采样因子，平衡计算量和精度
rx_signal_ds = rx_signal(:, 1:downsample_factor:end, 1:downsample_factor:end);
[~, n_chirps_ds, n_samples_ds] = size(rx_signal_ds);

% 重塑采样后的接收信号为矢量形式
y = reshape(rx_signal_ds, [], 1);

% 计算网格尺寸
n_R = length(R_grid);
n_theta = length(theta_grid);
n_phi = length(phi_grid);
n_atoms = n_R * n_theta * n_phi;

fprintf('字典参数: %d距离点 x %d方位角点 x %d俯仰角点 = %d个原子\n', ...
    n_R, n_theta, n_phi, n_atoms);

% ===== 内存高效的OMP实现 =====
% 不事先创建完整字典矩阵，而是按需计算

% 初始化
residual = y;
selected_params = [];
residual_norm = norm(residual);
initial_residual_norm = residual_norm;
iter = 0;

% 创建用于保存选中原子的矩阵
Phi_selected = [];

fprintf('开始OMP迭代，目标残差: %.2e\n', params.residual_thr * initial_residual_norm);
fprintf('信号向量y维度: %d\n', length(y));

% OMP迭代
while (iter < params.max_iterations) && (residual_norm > params.residual_thr * initial_residual_norm)
    iter = iter + 1;
    fprintf('OMP迭代 %d/%d，当前残差: %.2e\n', iter, params.max_iterations, residual_norm);
    
    % 找到与残差最相关的原子
    max_correlation = -1;
    best_atom = [];
    best_atom_params = [0, 0, 0]; % 初始化为默认值，避免空数组访问错误
    best_f_d = 0;
    
    % 预先计算距离、角度对应的权重
    R_weights = compute_distance_weights(R_grid, R_prior, sigma_R);
    theta_weights = compute_angle_weights(theta_grid, theta_prior, sigma_theta);
    phi_weights = compute_angle_weights(phi_grid, phi_prior, sigma_phi);
    
    % 进度显示参数
    progress_step = max(1, floor(n_R/5));
    
    % 对每个可能的参数组合计算相关性
    for i_R = 1:n_R
        R = R_grid(i_R);
        
        % 显示进度
        if mod(i_R, progress_step) == 0 || i_R == 1 || i_R == n_R
            fprintf('搜索距离 %.2f m (%d/%d - %.1f%%)\n', R, i_R, n_R, 100*i_R/n_R);
        end
        
        % 获取当前距离的权重
        w_R = R_weights(i_R);
        
        for i_theta = 1:n_theta
            theta_deg = theta_grid(i_theta);
            theta = theta_deg * pi/180;
            
            % 获取当前方位角的权重
            w_theta = theta_weights(i_theta);
            
            for i_phi = 1:n_phi
                phi_deg = phi_grid(i_phi);
                phi = phi_deg * pi/180;
                
                % 获取当前俯仰角的权重
                w_phi = phi_weights(i_phi);
                
                % 计算组合权重 (高斯先验) - 优化权重分配
                w = 0.7*w_R + 0.15*w_theta + 0.15*w_phi;  % 强调距离成分
                
                % 只考虑先验附近较强的权重，节省计算
                if w < 0.05
                    continue;
                end
                
                % 为每个速度假设生成字典原子并测试
                for v_idx = 1:length(v_r_options)
                    try
                        % 使用当前速度假设生成原子
                        current_f_d = f_d_options(v_idx);
                        atom = generate_atom_with_doppler(R, theta, phi, current_f_d, params, downsample_factor, residual);
                        
                        % 确保atom不为空
                        if isempty(atom) || any(isnan(atom))
                            continue;
                        end
                        
                        % 应用先验权重
                        atom = w * atom;
                        
                        % 计算与残差的相关性
                        correlation = abs(atom' * residual);
                        
                        % 更新最大相关性
                        if correlation > max_correlation
                            max_correlation = correlation;
                            best_atom = atom;
                            best_atom_params = [R, theta_deg, phi_deg];
                            best_f_d = current_f_d; % 记录最佳多普勒频移
                        end
                    catch e
                        % 出错时继续下一个速度假设
                        continue;
                    end
                end
            end
        end
    end
    
    % 打印最佳多普勒频移信息
    if max_correlation > 0
        v_r_est = best_f_d * params.c / (2 * params.fc);
        fprintf('  最佳速度估计: %.2f m/s (相关性: %.2e)\n', v_r_est, max_correlation);
    end
    
    % 检查是否找到有效原子
    if isempty(best_atom) || max_correlation <= 0
        fprintf('警告: 未找到有效原子，终止迭代\n');
        break;
    end
    
    % 将最佳原子参数添加到选定集合
    selected_params = [selected_params; best_atom_params];
    
    % 更新字典矩阵
    if isempty(Phi_selected)
        % 第一次迭代
        Phi_selected = best_atom;
    else
        % 后续迭代
        Phi_selected = [Phi_selected, best_atom];
    end
    
    % 最小二乘估计 - 使用QR分解提高数值稳定性
    try
        [Q, R_matrix] = qr(Phi_selected, 0);
        s = R_matrix \ (Q' * y);
    catch
        % 如果QR分解失败，使用伪逆
        s = pinv(Phi_selected) * y;
    end
    
    % 更新残差
    residual = y - Phi_selected * s;
    residual_norm = norm(residual);
    
    fprintf('已选择参数: R=%.2f m, θ=%.2f°, φ=%.2f°\n', ...
        best_atom_params(1), best_atom_params(2), best_atom_params(3));
    fprintf('剩余残差: %.2e (%.2f%%)\n', residual_norm, 100*residual_norm/initial_residual_norm);
end

% 检查是否成功找到任何参数
if isempty(selected_params)
    fprintf('未找到任何有效参数，返回先验估计\n');
    R_est = R_prior;
    theta_est = theta_prior;
    phi_est = phi_prior;
    return;
end

% 加权平均估计结果 (权重为稀疏系数的绝对值)
weights = abs(s);

% 检查权重是否有效
if sum(weights) > 0
    % 应用非线性变换，增强大权重的影响，抑制小权重的影响
    weights = weights.^1.5; % 使用幂函数增强权重差异
    
    % 移除异常值 - 如果某个权重远小于最大权重，将其设为0
    max_weight = max(weights);
    weights(weights < 0.05 * max_weight) = 0;
    
    % 重新归一化
    if sum(weights) > 0
        weights = weights / sum(weights);
    else
        % 如果所有权重都被移除，使用原始权重
        weights = abs(s);
        weights = weights / sum(weights);
    end
    
    % 计算加权平均
    R_est = sum(weights .* selected_params(:, 1));
    theta_est = sum(weights .* selected_params(:, 2));
    phi_est = sum(weights .* selected_params(:, 3));
    
    % 检查结果是否有效
    if isnan(R_est) || isnan(theta_est) || isnan(phi_est) || isinf(R_est) || isinf(theta_est) || isinf(phi_est)
        fprintf('警告: 加权平均结果包含NaN或Inf，使用最大权重对应的参数\n');
        [~, max_idx] = max(abs(s));
        R_est = selected_params(max_idx, 1);
        theta_est = selected_params(max_idx, 2);
        phi_est = selected_params(max_idx, 3);
    end
    
    % 检查估计结果是否在合理范围内
    if R_est <= 0 || R_est > 1000 || abs(theta_est) > 180 || abs(phi_est) > 90
        fprintf('警告: 估计结果超出合理范围，使用最后一次迭代的参数\n');
        R_est = selected_params(end, 1);
        theta_est = selected_params(end, 2);
        phi_est = selected_params(end, 3);
    end
else
    % 权重无效，使用最大系数对应的参数
    fprintf('警告: 无效权重，使用最后一次迭代的参数\n');
    R_est = selected_params(end, 1);
    theta_est = selected_params(end, 2);
    phi_est = selected_params(end, 3);
end
end

function weights = compute_distance_weights(grid, prior, sigma)
% 计算距离网格的权重
weights = exp(-(grid - prior).^2/(1.2*sigma^2));
end

function weights = compute_angle_weights(grid, prior, sigma)
% 计算角度网格的权重
weights = exp(-(grid - prior).^2/(2.5*sigma^2));
end

function grid = create_nonuniform_grid(min_val, max_val, center, step)
% 创建非均匀网格 - 中心附近更密集
if min_val >= max_val
    grid = center;
    return;
end

% 中心附近的密集网格
dense_min = max(min_val, center - 2*step);
dense_max = min(max_val, center + 2*step);
dense_grid = dense_min:(step*0.5):dense_max;

% 中心外围的稀疏网格
if min_val < dense_min
    sparse_grid_left = min_val:step:dense_min;
else
    sparse_grid_left = [];
end

if max_val > dense_max
    sparse_grid_right = (dense_max+step):step:max_val;
else
    sparse_grid_right = [];
end

% 组合网格并去除重复值
grid = unique([sparse_grid_left, dense_grid, sparse_grid_right]);
end