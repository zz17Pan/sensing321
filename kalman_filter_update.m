function kf = kalman_filter_update(kf, z, params)
% KALMAN_FILTER_UPDATE 更新卡尔曼滤波器状态
%   kf = KALMAN_FILTER_UPDATE(kf, z, params) 使用最新观测向量更新卡尔曼滤波器
%   z = [R; theta; phi] - 观测向量（距离，方位角，俯仰角，弧度制）

% 显示输入观测，确保单位和格式正确
fprintf('卡尔曼滤波器更新:\n');
fprintf('  观测: 距离=%.2fm, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
    z(1), z(2)*180/pi, z(3)*180/pi);

% 保存上一个状态以备调试
x_prev = kf.x;

% ===== 预测步骤 =====
% 1. 状态预测
x_pred = kf.F * kf.x;

% 2. 协方差预测
P_pred = kf.F * kf.P * kf.F' + kf.Q;

% 确保协方差矩阵是对称的
P_pred = (P_pred + P_pred')/2;

% ===== 计算当前状态下的预测位置和角度 =====
x = x_pred(1); % x位置
y = x_pred(4); % y位置
z_pos = x_pred(7); % z位置
vx = x_pred(2); % x速度
vy = x_pred(5); % y速度
vz = x_pred(8); % z速度

% 检查数值错误
if any(isnan([x, y, z_pos])) || any(isinf([x, y, z_pos]))
    fprintf('警告: 预测位置包含NaN或Inf，使用上一个有效状态\n');
    x = x_prev(1);
    y = x_prev(4);
    z_pos = x_prev(7);
    vx = x_prev(2);
    vy = x_prev(5);
    vz = x_prev(8);
    % 更新预测状态向量
    x_pred(1) = x;
    x_pred(4) = y;
    x_pred(7) = z_pos;
    x_pred(2) = vx;
    x_pred(5) = vy;
    x_pred(8) = vz;
end

% 计算预测的距离 - 确保不为零
r_pred = sqrt(x^2 + y^2 + z_pos^2 + 1e-10); % 添加小常数避免除零

% 计算预测的方位角，标准定义：逆时针，相对x轴
theta_pred = atan2(y, x);

% 计算预测的俯仰角，确保分母不为零
phi_pred = atan2(z_pos, sqrt(x^2 + y^2 + 1e-10));

% 报告预测值
fprintf('  状态预测: 位置=[%.2f, %.2f, %.2f]m, 速度=[%.2f, %.2f, %.2f]m/s\n', ...
    x, y, z_pos, vx, vy, vz);
fprintf('  预测观测: 距离=%.2fm, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
    r_pred, theta_pred*180/pi, phi_pred*180/pi);

% ===== 计算雅可比矩阵 =====
% 观测方程对状态的偏导数 - 使用数值稳定的实现
% 为避免除以零，添加小常数
epsilon = 1e-10;

% 距离偏导数
dr_dx = x / (r_pred + epsilon);
dr_dy = y / (r_pred + epsilon);
dr_dz = z_pos / (r_pred + epsilon);

% 方位角偏导数，确保分母不为零
xy_sq = x^2 + y^2 + epsilon;
dtheta_dx = -y / xy_sq;
dtheta_dy = x / xy_sq;
dtheta_dz = 0;

% 俯仰角偏导数，确保所有分母都不为零
r_xy = sqrt(x^2 + y^2 + epsilon);
dphi_dx = -x * z_pos / ((r_pred^2) * r_xy + epsilon);
dphi_dy = -y * z_pos / ((r_pred^2) * r_xy + epsilon);
dphi_dz = r_xy / (r_pred^2 + epsilon);

% 构建雅可比矩阵 (3x9)
H = zeros(3, 9);
H(1, 1) = dr_dx;     H(1, 4) = dr_dy;     H(1, 7) = dr_dz;
H(2, 1) = dtheta_dx; H(2, 4) = dtheta_dy; H(2, 7) = dtheta_dz;
H(3, 1) = dphi_dx;   H(3, 4) = dphi_dy;   H(3, 7) = dphi_dz;

% 检查雅可比矩阵的数值问题
if any(isnan(H(:))) || any(isinf(H(:)))
    fprintf('警告: 雅可比矩阵包含NaN或Inf，使用单位矩阵替代\n');
    % 构造简化的雅可比矩阵，确保数值稳定
    H = zeros(3, 9);
    H(1, 1) = 1/r_pred; H(1, 4) = 1/r_pred; H(1, 7) = 1/r_pred;  % 距离对位置的偏导数估计
    H(2, 1) = 1/xy_sq;  H(2, 4) = 1/xy_sq;  % 方位角对x,y的偏导数估计
    H(3, 7) = 1/r_pred; % 俯仰角对z的偏导数估计
end

% ===== 更新步骤 =====
% 0. 构建观测向量和预测观测
z_pred = [r_pred; theta_pred; phi_pred];

% 1. 计算创新向量（测量残差）
innovation = z - z_pred;

% 1.1 处理角度环绕（确保角度差在-pi和pi之间）
innovation(2) = wrapToPi(innovation(2));
innovation(3) = wrapToPi(innovation(3));

% 1.2 检查创新向量是否有异常值，限制创新大小
% 优化匀速运动场景下的创新限制
max_innovation = [2.0; 6*pi/180; 6*pi/180]; % 减小最大允许创新量：距离2.0m，角度6度

% 应用自适应创新限制 - 根据当前状态协方差动态调整
% 对于不确定性高的状态，允许更大的创新；对于不确定性低的状态，限制创新
pos_uncertainty = sqrt(diag(P_pred([1,4,7], [1,4,7]))); % 位置状态的不确定性
vel_uncertainty = sqrt(diag(P_pred([2,5,8], [2,5,8]))); % 速度状态的不确定性

% 改进的自适应因子计算，使其更平滑且考虑速度不确定性
% 当速度不确定性大时，适当增大创新限制
pos_factor = min(1.0, max(0.3, 0.35 * norm(pos_uncertainty) / 8.0));
vel_factor = min(0.3, max(0.1, 0.2 * norm(vel_uncertainty) / 5.0));
scale_factor = pos_factor + vel_factor;

% 动态调整最大创新量 - 使用更保守的基础值
max_innovation = max_innovation * (0.55 + scale_factor);

% 应用创新限制
innovation = min(max(innovation, -max_innovation), max_innovation);

fprintf('  创新向量: 距离=%.2fm, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
    innovation(1), innovation(2)*180/pi, innovation(3)*180/pi);

% 2. 计算创新协方差
S = H * P_pred * H' + kf.R;

% 2.1 确保S是对称正定的
S = (S + S')/2;
[V, D] = eig(S);
D = diag(D);
D = max(D, 1e-6); % 确保最小特征值大于零
S = V * diag(D) * V';

% 3. 计算卡尔曼增益
K = P_pred * H' / S; % 直接使用矩阵除法，更加稳定

% 3.1 应用自适应增益调整 - 根据创新大小动态调整增益
% 当创新过大时，适当减小增益，提高稳定性
% 使用加权归一化，给距离和角度不同的权重
innovation_weights = [0.8; 1.0*pi/180; 1.0*pi/180]; % 降低距离权重，增加角度权重的相对影响
innovation_norm = norm(innovation ./ innovation_weights); % 加权归一化创新向量

% 计算速度变化率 - 用于判断运动状态
vel_change = norm([kf.x(2)-x_prev(2); kf.x(5)-x_prev(5); kf.x(8)-x_prev(8)]);
vel_magnitude = norm([kf.x(2); kf.x(5); kf.x(8)]);
vel_change_ratio = 0;
if vel_magnitude > 0.1
    vel_change_ratio = vel_change / vel_magnitude;
end

% 使用更平滑的增益调整函数，考虑运动状态
if innovation_norm > 3.0
    % 使用平滑的非线性函数降低增益，避免突变
    % 当速度变化率小时(匀速运动)，使用更小的增益
    base_scale = max(0.3, min(0.9, 3.0 / innovation_norm));
    if vel_change_ratio < 0.1 % 匀速运动
        gain_scale = base_scale * 0.85; % 进一步降低增益
    else
        gain_scale = base_scale;
    end
    fprintf('  应用增益缩放: %.2f (创新范数: %.2f, 速度变化率: %.2f)\n', gain_scale, innovation_norm, vel_change_ratio);
    K = K * gain_scale;
elseif innovation_norm < 0.6
    % 当创新很小时，可以适当增大增益以加快收敛
    % 但在匀速运动时保持较小增益
    if vel_change_ratio < 0.1 % 匀速运动
        gain_scale = min(1.05, 1.0 + 0.05*(0.6 - innovation_norm));
    else
        gain_scale = min(1.1, 1.0 + 0.1*(0.6 - innovation_norm));
    end
    fprintf('  应用增益增强: %.2f (创新范数: %.2f, 速度变化率: %.2f)\n', gain_scale, innovation_norm, vel_change_ratio);
    K = K * gain_scale;
end

% 检查卡尔曼增益是否包含NaN或Inf
if any(isnan(K(:))) || any(isinf(K(:)))
    fprintf('警告: 卡尔曼增益包含NaN或Inf，使用小常数值替代\n');
    K = 0.1 * ones(size(K)); % 使用小常数作为安全值
end

% 4. 更新状态
kf.x = x_pred + K * innovation;

% 5. 更新状态协方差矩阵 (使用Joseph稳定形式)
I = eye(size(kf.P));
kf.P = (I - K * H) * P_pred * (I - K * H)' + K * kf.R * K';

% 5.0.1 应用协方差下限 - 防止协方差过小导致滤波器过度自信
% 增大最小方差值，防止滤波器过度自信
min_position_var = 0.15^2; % 位置最小方差 - 增大以防止过度自信
min_velocity_var = 0.6^2;  % 速度最小方差 - 增大以适应速度变化
min_accel_var = 0.8^2;     % 加速度最小方差 - 减小以适应匀速运动

% 计算当前速度大小，用于自适应调整协方差下限
vel_magnitude = norm([kf.x(2); kf.x(5); kf.x(8)]);

% 自适应调整协方差下限 - 速度越大，位置不确定性越大
vel_factor = min(2.0, max(1.0, vel_magnitude / 10.0));
adaptive_pos_var = min_position_var * vel_factor;

% 确保位置状态的协方差不会过小
kf.P(1,1) = max(kf.P(1,1), adaptive_pos_var);
kf.P(4,4) = max(kf.P(4,4), adaptive_pos_var);
kf.P(7,7) = max(kf.P(7,7), adaptive_pos_var);

% 确保速度状态的协方差不会过小
kf.P(2,2) = max(kf.P(2,2), min_velocity_var);
kf.P(5,5) = max(kf.P(5,5), min_velocity_var);
kf.P(8,8) = max(kf.P(8,8), min_velocity_var);

% 确保加速度状态的协方差不会过小
kf.P(3,3) = max(kf.P(3,3), min_accel_var);
kf.P(6,6) = max(kf.P(6,6), min_accel_var);
kf.P(9,9) = max(kf.P(9,9), min_accel_var);

% 打印协方差下限信息
fprintf('  协方差下限: 位置=%.4f, 速度=%.4f, 加速度=%.4f\n', adaptive_pos_var, min_velocity_var, min_accel_var);

% 5.1 确保P是对称正定的
kf.P = (kf.P + kf.P')/2;
[V, D] = eig(kf.P);
D = diag(D);
D = max(D, 1e-6); % 确保最小特征值大于零
kf.P = V * diag(D) * V';

% ===== 报告更新后的状态 =====
x = kf.x(1); 
y = kf.x(4);
z_pos = kf.x(7);
vx = kf.x(2);
vy = kf.x(5);
vz = kf.x(8);
ax = kf.x(3);
ay = kf.x(6);
az = kf.x(9);

% 检查更新后的状态是否合理
if any(isnan([x, y, z_pos])) || any(isinf([x, y, z_pos]))
    fprintf('警告: 更新后状态包含NaN或Inf，使用预测状态替代\n');
    kf.x = x_pred;
    x = kf.x(1); 
    y = kf.x(4);
    z_pos = kf.x(7);
    vx = kf.x(2);
    vy = kf.x(5);
    vz = kf.x(8);
    ax = kf.x(3);
    ay = kf.x(6);
    az = kf.x(9);
end

fprintf('  更新后状态: 位置=[%.2f, %.2f, %.2f]m, 速度=[%.2f, %.2f, %.2f]m/s\n', ...
    x, y, z_pos, vx, vy, vz);
fprintf('  加速度: [%.2f, %.2f, %.2f]m/s²\n', ax, ay, az);
fprintf('  速度比较: 更新前=[%.2f, %.2f, %.2f]m/s, 更新后=[%.2f, %.2f, %.2f]m/s\n', ...
    x_prev(2), x_prev(5), x_prev(8), vx, vy, vz);

% 计算更新后的球坐标，确保不会出现除零错误
r_est = sqrt(x^2 + y^2 + z_pos^2 + 1e-10);
theta_est = atan2(y, x);
phi_est = atan2(z_pos, sqrt(x^2 + y^2 + 1e-10));

fprintf('  对应球坐标: 距离=%.2fm, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
    r_est, theta_est*180/pi, phi_est*180/pi);

end

function angle = wrapToPi(angle)
% 将角度包装到-pi到pi的范围内
    angle = mod(angle + pi, 2*pi) - pi;
end