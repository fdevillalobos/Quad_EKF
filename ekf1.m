function [X, Z] = ekf1(sensor, vic, K, Tb_to_c, Real)
% EKF1 Extended Kalman Filter with Vicon velocity as inputs
%
% INPUTS:
%   sensor - struct stored in provided dataset, fields include
%          - is_ready: logical, indicates whether sensor data is valid
%          - t: sensor timestamp
%          - rpy, omg, acc: imu readings, you should not use these in this phase
%          - img: uint8, 240x376 grayscale image
%          - id: 1xn ids of detected tags
%          - p0, p1, p2, p3, p4: 2xn pixel position of center and
%                                four corners of detected tags
%            Y
%            ^ P3 == P2
%            | || P0 ||
%            | P4 == P1
%            o---------> X
%   vic    - struct for storing vicon linear velocity in world frame and
%            angular velocity in body frame, fields include
%          - t: vicon timestamp
%          - vel = [vx; vy; vz; wx; wy; wz]
%   varargin - any variables you wish to pass into the function, could be
%              a data structure to represent the map or camera parameters,
%              your decision. But for the purpose of testing, since we don't
%              know what inputs you will use, you have to specify them in
%              init_script by doing
%              ekf1_handle = ...
%                  @(sensor, vic) ekf1(sensor, vic, your input arguments);
%
% OUTPUTS:
% X - nx1 state of the quadrotor, n should be greater or equal to 6
%     the state should be in the following order
%     [x; y; z; roll; pitch; yaw; other states you use]
%     we will only take the first 6 rows of X
% OPTIONAL OUTPUTS:
% Z - mx1 measurement of your pose estimator, m shoulb be greater or equal to 6
%     the measurement should be in the following order
%     [x; y; z; roll; pitch; yaw; other measurement you use]
%     note that this output is optional, it's here in case you want to log your
%     measurement for your final report for part 2
persistent Sigma_old mu_old W R Q t_old Z_old first_loop
estimate_pose_handle = @(sensor) estimate_pose(sensor, K, Tb_to_c, Real);

% Get z out of our pose estimation
[pos, eul] = estimate_pose_handle(sensor);
z = [pos; eul];

if(isempty(Sigma_old))
    delta_t   = 0.01;
    mu_old    = zeros(6,1);
    Sigma_old = eye(6) * 0.1;
    W         = eye(6);
    Q         = eye(6) * 0.5;
    %R = covariance of the data.
    R = eye(6) * 0.0001;
    Z_old     = zeros(6,1);
    first_flag = 1;
else
    delta_t = vic.t - t_old;
    first_flag = 0;
end

t_old = vic.t;

%%
% vicon_data = [x y z roll pitch yaw vx vy vz wx wy wz]'
H  = eye(6);
G  = G_function(mu_old, vic.vel, zeros(6,1), delta_t);
L  = L_function(mu_old, vic.vel, zeros(6,1), delta_t);
mu_bar = g_func(mu_old, vic.vel, zeros(6,1), delta_t);
Sigma_bar = G * Sigma_old * G' + L * Q * L';

if(first_loop == 0)
    mu_old = [pos; eul];
    Z = z;
    Sigma_old = Sigma_bar;
    first_loop = 1;
    X = [];
else
    if(~isempty(sensor) && ~isempty(sensor.id))
        k = Sigma_bar * H / (H * Sigma_bar * H' + W * R * W');
        mu_old = mu_bar + k * (z - H * mu_bar);
        Sigma_old = (eye(6) - k * H) * Sigma_bar;
        Z = z;
        Z_old = Z;
        X = mu_old;
    else
        X = mu_bar;
        Z = z;
        mu_old = mu_bar;
        Sigma_old = Sigma_bar;
    end
       
end

end
