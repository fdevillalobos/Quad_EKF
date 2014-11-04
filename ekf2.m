function [X, Z] = ekf2(sensor, K, Tb_to_c, Real)
% EKF2 Extended Kalman Filter with IMU as inputs
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
%   varargin - any variables you wish to pass into the function, could be
%              a data structure to represent the map or camera parameters,
%              your decision. But for the purpose of testing, since we don't
%              know what inputs you will use, you have to specify them in
%              init_script by doing
%              ekf1_handle = ...
%                  @(sensor) ekf2(sensor, your input arguments);
%
% OUTPUTS:
% X - nx1 state of the quadrotor, n should be greater or equal to 9
%     the state should be in the following order
%     [x; y; z; vx; vy; vz; roll; pitch; yaw; other states you use]
%     we will only take the first 9 rows of X
% OPTIONAL OUTPUTS:
% Z - mx1 measurement of your pose estimator, m shoulb be greater or equal to 6
%     the measurement should be in the following order
%     [x; y; z; roll; pitch; yaw; other measurement you use]
%     note that this output is optional, it's here in case you want to log your
%     measurement for your final report for part 2

persistent Sigma_old mu_old W R Q t_old Z_old
estimate_pose_handle = @(sensor) estimate_pose(sensor, K, Tb_to_c, Real);

if(isempty(Sigma_old))
    [pos, eul] = estimate_pose_handle(sensor);
    delta_t   = 0.01;
    vel_i = [0.0;  0.0; 0.0];
    biasa = [0.1321; 0.0216; 9.724]; %-0.086];
    biasg = [-1.2e-05; -2.6e-05; -0.0011];
    mu_old    =  [pos; eul;  vel_i; biasa; biasg];
    Sigma_old =  eye(15) * 0.1;
    W         = -eye(6);
    Q         =  eye(12) * 1;
    
    %R = covariance of the data.
    R1 = 0.001 * [0.0048  0.0089 -0.0015;
                  0.0089  0.9057  0.0541;
                  -0.0015 0.0541 0.0042];
              
    R2 = 0.0001* [ 0.7639  0.1845 -0.5129;
                   0.1845  0.4423 -0.3234;
                  -0.5129 -0.3234 0.6025];
    R = [R1, zeros(3);
        zeros(3), R2];
%     R = 0.001 * eye(6);
    
    Z_old     = zeros(6,1);
else
    delta_t = sensor.t - t_old;
end
    
t_old = sensor.t;
    
    %%
    % vicon_data = [x y z roll pitch yaw vx vy vz wx wy wz]'
    H  = [eye(6) zeros(6,9)];
    G  = G2_function(mu_old, sensor.omg, sensor.acc, zeros(6,1), zeros(6,1), delta_t);
    L  = L2_function(mu_old, sensor.omg, sensor.acc, zeros(6,1), zeros(6,1), delta_t);
    mu_bar = g2_func(mu_old, sensor.omg, sensor.acc, zeros(6,1), zeros(6,1), delta_t, 9.81);
    Sigma_bar = G * Sigma_old * G' + L * Q * L';
   
    if(~isempty(sensor) && ~isempty(sensor.id))
        % % % Get z out of our pose estimation
        [pos, eul] = estimate_pose_handle(sensor);
        z = [pos; eul];
        k = Sigma_bar * H' / (H * Sigma_bar * H' + W * R * W');
        mu_old = mu_bar + k * (z - H * mu_bar);
        Sigma_old = (eye(15) - k * H) * Sigma_bar;                          
        Z = z;
        Z_old = Z;
        X = mu_old(1:9);
        Aux = X(4:6);
        X(4:6) = X(7:9);
        X(7:9) = Aux;
    else
        Z = Z_old;
        mu_old = mu_bar;
        X = [];
        Sigma_old = Sigma_bar;
    end  
    
end
