% Add additional inputs after the given ones if you want to
% Example:
% your_input = 1;
% ekf_handle1 = @(sensor, vic) ekf1(sensor, vic, your_input);
% ekf_handle2 = @(sensor) ekf2(sensor, your_input);
%
% We will only call ekf_handle in the test function.
% Note that this will only create a function handle, but not run the function

%% Initialize Handlers and data
% clear all; close all; clc;

%% Calculate the World coordinates of the April Tags.
% Camera Matrix (zero-indexed):
K = [314.1779       0       199.4848; ...
          0      314.2218   113.7838; ...
          0         0          1];

% Camera-IMU Calibration (see attached images for details):
Tb_to_c = [-0.04; 0.0; -0.03];

% Calculate the real world coord of the data.
rows = 12;      cols = 9;
x_dim = 0.152;
y_dim = 0.152;
y_div = [0.152 0.152 0.178 0.152 0.152 0.178 0.152 0.152];
x_div = 0.152 * ones(11,1)';

Real = world(x_dim, y_dim, x_div, y_div, rows, cols);

% ekf1_jacobians;
%% Initialize Handlers and data
ekf1_handle = @(sensor, vic) ekf1(sensor, vic, K, Tb_to_c, Real);
ekf2_handle = @(sensor) ekf2(sensor, K, Tb_to_c, Real);

%%
%test_script_2;











