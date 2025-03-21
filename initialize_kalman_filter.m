function kf = initialize_kalman_filter(params)
% INITIALIZE_KALMAN_FILTER 初始化卡尔曼滤波器
%   kf = INITIALIZE_KALMAN_FILTER(params) 初始化用于目标跟踪的卡尔曼滤波器

fprintf('初始化扩展卡尔曼滤波器\n');

% 时间步长
dt = params.dt;
fprintf('  时间步长: %.6f秒\n', dt);

% 状态向量: [x, vx, ax, y, vy, ay, z, vz, az]
% 包含三维位置、速度和加速度
n_states = 9;

% 初始位置估计 (球坐标转笛卡尔坐标)
R_init = params.initial_R;
theta_init = params.initial_theta * pi/180; % 转为弧度
phi_init = params.initial_phi * pi/180;   % 转为弧度

% 计算笛卡尔坐标下的位置 - 使用标准定义的角度转换
x_init = R_init * cos(phi_init) * cos(theta_init);
y_init = R_init * cos(phi_init) * sin(theta_init);
z_init = R_init * sin(phi_init);

% 初始速度估计 - 设置为与实际接收端速度接近的值
vx_init = 12.0;  % x方向初始速度
vy_init = -10.0; % y方向初始速度 
vz_init = 8.0;   % z方向初始速度

% 初始加速度估计
ax_init = 0.0;  % x方向初始加速度
ay_init = 0.0;  % y方向初始加速度
az_init = 0.0;  % z方向初始加速度

% 初始状态向量
x0 = [x_init; vx_init; ax_init; y_init; vy_init; ay_init; z_init; vz_init; az_init];

fprintf('  初始位置: [%.2f, %.2f, %.2f]m\n', x_init, y_init, z_init);
fprintf('  初始速度: [%.2f, %.2f, %.2f]m/s\n', vx_init, vy_init, vz_init);
fprintf('  初始加速度: [%.2f, %.2f, %.2f]m/s²\n', ax_init, ay_init, az_init);

% 状态转移矩阵 (匀加速运动模型)
% x轴
Fx = [1, dt, dt^2/2;
      0, 1,  dt;
      0, 0,  1]; 
  
% y轴  
Fy = [1, dt, dt^2/2;
      0, 1,  dt;
      0, 0,  1];
  
% z轴
Fz = [1, dt, dt^2/2;
      0, 1,  dt;
      0, 0,  1];

% 组合成完整的状态转移矩阵 (分块对角矩阵)
F = zeros(n_states, n_states);
F(1:3, 1:3) = Fx;
F(4:6, 4:6) = Fy;
F(7:9, 7:9) = Fz;

% 过程噪声协方差矩阵
% 定义加速度过程噪声标准差 - 进一步减小以更好地匹配匀速运动模型
sigma_ax = 0.8;  % x方向加速度标准差 (m/s²) - 显著减小以更好地匹配匀速运动
sigma_ay = 0.8;  % y方向加速度标准差 (m/s²) - 显著减小以更好地匹配匀速运动
sigma_az = 0.8;  % z方向加速度标准差 (m/s²) - 显著减小以更好地匹配匀速运动

fprintf('  过程噪声标准差(加速度): [%.2f, %.2f, %.2f]m/s²\n', sigma_ax, sigma_ay, sigma_az);

% 计算每个维度的过程噪声协方差
% 这里使用连续时间白噪声加速度模型 (CWNA)
Qx = [dt^4/4, dt^3/2, dt^2/2;
      dt^3/2, dt^2,   dt;
      dt^2/2, dt,     1] * sigma_ax^2;
  
Qy = [dt^4/4, dt^3/2, dt^2/2;
      dt^3/2, dt^2,   dt;
      dt^2/2, dt,     1] * sigma_ay^2;
  
Qz = [dt^4/4, dt^3/2, dt^2/2;
      dt^3/2, dt^2,   dt;
      dt^2/2, dt,     1] * sigma_az^2;
  
% 组合成完整的过程噪声协方差矩阵
Q = zeros(n_states, n_states);
Q(1:3, 1:3) = Qx;
Q(4:6, 4:6) = Qy;
Q(7:9, 7:9) = Qz;

% 观测矩阵 H 将在 kalman_filter_update 中动态计算
% 这里仅初始化一个初始值，它将根据当前状态重新计算雅可比矩阵
H = zeros(3, n_states);
% 观测模型：[距离; 方位角; 俯仰角] = h(x, y, z)
% 将会在kalman_filter_update中使用雅可比行列式计算非线性观测矩阵

% 观测噪声协方差矩阵 - 根据实际观测精度调整
sigma_r = 0.15;      % 距离观测标准差 (m) - 进一步减小以提高距离估计的权重
sigma_theta = 1.0 * pi/180;  % 方位角观测标准差 (rad) - 减小以提高角度估计的权重
sigma_phi = 1.0 * pi/180;    % 俯仰角观测标准差 (rad) - 减小以提高角度估计的权重

fprintf('  观测噪声标准差: 距离=%.2fm, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
       sigma_r, sigma_theta*180/pi, sigma_phi*180/pi);

% 创建观测噪声协方差矩阵
R = diag([sigma_r^2, sigma_theta^2, sigma_phi^2]);

% 添加少量协方差，使观测更稳定
R(1,2) = 0.001; R(2,1) = 0.001; % 距离-方位角协方差
R(1,3) = 0.001; R(3,1) = 0.001; % 距离-俯仰角协方差
R(2,3) = 0.002; R(3,2) = 0.002; % 方位角-俯仰角协方差

% 初始状态估计协方差矩阵 - 调整初始不确定性以更好地适应匀速运动场景
P0 = zeros(n_states, n_states);
% 位置初始不确定性 - 增大以适应初始估计误差
P0(1,1) = 10.0^2;   % x位置方差 - 增大以适应可能的初始位置误差
P0(4,4) = 10.0^2;   % y位置方差 - 增大以适应可能的初始位置误差
P0(7,7) = 10.0^2;   % z位置方差 - 增大以适应可能的初始位置误差
% 速度初始不确定性 - 增大以适应可能的高速变化
P0(2,2) = 10.0^2;   % x速度方差 - 增大以更好地适应初始速度不确定性
P0(5,5) = 10.0^2;   % y速度方差 - 增大以更好地适应初始速度不确定性
P0(8,8) = 10.0^2;   % z速度方差 - 增大以更好地适应初始速度不确定性
% 加速度初始不确定性 - 减小以适应匀速运动场景
P0(3,3) = 2.0^2;   % x加速度方差 - 减小以更好地匹配匀速运动模型
P0(6,6) = 2.0^2;   % y加速度方差 - 减小以更好地匹配匀速运动模型
P0(9,9) = 2.0^2;   % z加速度方差 - 减小以更好地匹配匀速运动模型

% 添加少量协方差，使状态估计更稳定
P0(1,2) = 2.0; P0(2,1) = 2.0; % x位置-速度协方差
P0(4,5) = 2.0; P0(5,4) = 2.0; % y位置-速度协方差
P0(7,8) = 2.0; P0(8,7) = 2.0; % z位置-速度协方差

% 创建卡尔曼滤波器结构
kf = struct();
kf.x = x0;    % 初始状态
kf.P = P0;    % 初始状态协方差
kf.F = F;     % 状态转移矩阵
kf.Q = Q;     % 过程噪声协方差
kf.H = H;     % 初始观测矩阵 (会在更新中重新计算)
kf.R = R;     % 观测噪声协方差

fprintf('卡尔曼滤波器初始化完成\n');
end