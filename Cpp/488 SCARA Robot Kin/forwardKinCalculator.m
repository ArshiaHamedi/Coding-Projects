% MATLAB script for RRPR robot forward kinematics
clear all;
clc;

% Prompt user for input parameters
theta1 = input('Enter theta1 (in degrees): ');
theta2 = input('Enter theta2 (in degrees): ');
d3 = input('Enter d3 (prismatic joint displacement): ');
theta4 = input('Enter theta4 (in degrees): ');

% Convert angles from degrees to radians
theta1 = deg2rad(theta1);
theta2 = deg2rad(theta2);
theta4 = deg2rad(theta4);

% Define DH parameters
% For RRPR robot: 2 revolute joints, 1 prismatic joint, 1 revolute joint
% Assuming unit lengths for links (L1 = 1, L2 = 1) for simplicity
% You can modify these lengths as needed
L1 = 1;  % Length of first link
L2 = 1;  % Length of second link

% DH Parameter Table
% [theta, d, a, alpha]
dh_params = [
    theta1     0    L1    0;    % Joint 1 (Revolute)
    theta2     0    L2    0;    % Joint 2 (Revolute)
    0         d3    0     0;    % Joint 3 (Prismatic)
    theta4     0    0     0     % Joint 4 (Revolute)
];

% Initialize the transformation matrix
T_final = eye(4);

% Calculate the homogeneous transformation matrix for each joint
for i = 1:4
    theta = dh_params(i,1);
    d = dh_params(i,2);
    a = dh_params(i,3);
    alpha = dh_params(i,4);
    
    % Individual transformation matrix (Denavit-Hartenberg)
    T = [
        cos(theta)  -sin(theta)*cos(alpha)  sin(theta)*sin(alpha)   a*cos(theta);
        sin(theta)   cos(theta)*cos(alpha) -cos(theta)*sin(alpha)   a*sin(theta);
        0            sin(alpha)             cos(alpha)              d;
        0            0                      0                       1
    ];
    
    % Multiply to get the final transformation
    T_final = T_final * T;
end

% Display results
fprintf('\nForward Kinematics Transformation Matrix (0 to 4):\n');
disp(T_final);

% Extract and display position and orientation
px = T_final(1,4);
py = T_final(2,4);
pz = T_final(3,4);

fprintf('End-effector position:\n');
fprintf('px = %.4f\n', px);
fprintf('py = %.4f\n', py);
fprintf('pz = %.4f\n', pz);

% For a planar robot, we can extract the total orientation angle
total_angle = rad2deg(atan2(T_final(2,1), T_final(1,1)));
fprintf('End-effector orientation (degrees): %.4f\n', total_angle);