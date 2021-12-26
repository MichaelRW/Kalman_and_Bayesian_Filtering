%This file is for Problem 7.6, Part (c)

clear all
 
randn('seed',12);           % Choose a randomization seed
xstate=[100; 0.3];          % State vector
Pcov=[45 0.2; 0.2 0.001];   % Error covariance matrix
R=[1 0; 0 1];               % Measurement noise matrix
 
% Linearize the nonlinear measurement model by computing H matrix of
% partial derivatives
H=[cos(xstate(2,1)) -xstate(1,1)*sin(xstate(2,1)); sin(xstate(2,1)) xstate(1,1)*cos(xstate(2,1))];
rescov=H*Pcov*H'+R;     % Standard Kalman filter measurement residual covariance
Kgain=Pcov*H'/rescov;   % Standard Kalman filter gain
disp(rescov);
disp(Kgain);

