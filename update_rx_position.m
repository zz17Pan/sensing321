function rx_array = update_rx_position(rx_array, params, frame_idx)
% UPDATE_RX_POSITION 更新接收端天线阵列的位置
%   rx_array = UPDATE_RX_POSITION(rx_array, params, frame_idx) 
%   根据运动模型更新接收端位置和角度

% 时间增量
dt = params.dt;

fprintf('\n帧 %d: 位置更新\n', frame_idx);
fprintf('  更新前位置: [%.2f, %.2f, %.2f]m\n', rx_array.position(1), rx_array.position(2), rx_array.position(3));
fprintf('  速度矢量: [%.2f, %.2f, %.2f]m/s, 时间步长: %.6fs\n', ...
    rx_array.velocity(1), rx_array.velocity(2), rx_array.velocity(3), dt);

% 保存旧角度用于计算变化量
theta_old = rx_array.theta;
phi_old = rx_array.phi;

% 更新接收端整体位置 (三维匀速直线运动)
position_old = rx_array.position;
position_delta = rx_array.velocity * dt;
rx_array.position = position_old + position_delta;

fprintf('  位置变化: [%.2f, %.2f, %.2f]m\n', position_delta(1), position_delta(2), position_delta(3));
fprintf('  更新后位置: [%.2f, %.2f, %.2f]m\n', rx_array.position(1), rx_array.position(2), rx_array.position(3));

% 根据新位置计算角度 - 使用标准定义，确保与MUSIC角度定义一致
% 方位角(theta) - x-y平面内，相对于x轴正方向，逆时针为正
% 与MUSIC算法中的定义一致: u = cos(phi) * sin(theta)
rx_array.theta = atan2(rx_array.position(2), rx_array.position(1)) * 180/pi;

% 俯仰角(phi) - 与x-y平面的夹角，向上为正
% 与MUSIC算法中的定义一致: v = sin(phi)
rx_array.phi = atan2(rx_array.position(3), ...
                     sqrt(rx_array.position(1)^2 + rx_array.position(2)^2)) * 180/pi;

fprintf('  方位角: %.2f° → %.2f° (变化: %.2f°)\n', ...
    theta_old, rx_array.theta, rx_array.theta - theta_old);
fprintf('  俯仰角: %.2f° → %.2f° (变化: %.2f°)\n', ...
    phi_old, rx_array.phi, rx_array.phi - phi_old);

% 将角度转换为弧度，用于后续计算
theta_rad = rx_array.theta * pi/180;
phi_rad = rx_array.phi * pi/180;

% 更新子阵位置，使用与MUSIC一致的方向余弦
% 方向余弦定义：
% u = cos(phi) * sin(theta)
% v = sin(phi)
% w = cos(phi) * cos(theta)

% 生成旋转矩阵
% 先绕z轴旋转theta角度
Rz = [cos(theta_rad), -sin(theta_rad), 0;
      sin(theta_rad), cos(theta_rad),  0;
      0,              0,               1];

% 再绕x轴旋转phi角度
Rx = [1, 0,           0;
      0, cos(phi_rad), -sin(phi_rad);
      0, sin(phi_rad), cos(phi_rad)];

% 旋转矩阵
R = Rx * Rz;

% 更新子阵位置
for k = 1:length(rx_array.subarrays)
    % 子阵中心相对于阵列中心的位置偏移
    offset = [(k - 2.5) * params.d_sub, 0, 0];
    
    % 应用旋转
    rotated_offset = (R * offset')';
    
    % 更新子阵中心位置
    rx_array.subarrays{k}.center = rx_array.position + rotated_offset;
    
    % 更新子阵内天线位置
    for nx = 1:4
        for nz = 1:4
            % 计算天线相对于子阵中心的偏移
            antenna_offset = [(nx - 2.5) * params.d, 0, (nz - 2.5) * params.d];
            
            % 应用旋转
            rotated_antenna_offset = (R * antenna_offset')';
            
            % 更新天线位置
            rx_array.subarrays{k}.antenna_positions(nx, nz, :) = ...
                rx_array.subarrays{k}.center + rotated_antenna_offset;
        end
    end
end

end 