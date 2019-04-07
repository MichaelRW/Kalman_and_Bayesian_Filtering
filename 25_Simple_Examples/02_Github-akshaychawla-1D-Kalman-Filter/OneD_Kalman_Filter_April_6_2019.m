

%% Environment

close all; clear; clc;
set(0, 'DefaultFigureWindowStyle', 'docked');



%% System Variables

X = 0;
dt = 1;
u = 2; %speed/control = 2m/s
rng(0); n=randn(); v=randn();



%% Kalman Filter Variables

x = 50;      %state vector
A = 1;    %state transition matrix
B = 1;    %control input matrix 
P = 100;    %std_dev*std_dev = 10*10 

Q = 100;  % Process noise variance.
R = 9;  % Measurement noise variance.

H = 1;



%% Initialize Processing Variables

simulationSize=20;

XX = zeros(1, simulationSize);
tt = zeros(1, simulationSize);
xx_predicted = zeros(1, simulationSize);
xx = zeros(1, simulationSize);
PP = zeros(1, simulationSize);
yy = zeros(1, simulationSize);



%% Simulating the System

for t=1:1:simulationSize
    
    n = sqrt(Q) * randn();  % State prediction noise (zero mean Gaussian distribution).
        X = X + u*dt + n;  % Predict next state.
        
    v = sqrt(R) * randn();  % State measurement noise (zero mean Gaussian distribution).
        y = H*X + v;  % Measurement of state.
    
   % Prediction Step
    x_predicted = A*x + B*u;  % Predicting the new state (mean).
        P = A*P*A.' + Q;  % Update the uncertainty of the new state.
    
   % Correction Step
    e = H*x_predicted;  % Expectation:  predicted measurement from the o/p.
        E = H*P*H';  % Covariance of expectation.
    
    z = y - e;  % Innovation:  difference between the expected state and the measurement.
    Z = R + E;  % Covariance:  sum of uncertainties of expectation and the measurement.
    K = P*H' * Z^-1;  % Update the Kalman gain.
    
    x = x_predicted + K*z;  % Final corrected state.
        P = P - K * H* P;  % Update the uncertainity in corrected state.
    
   
    xx_predicted(t) = x_predicted;
    xx(t) = x; PP(t) = P; XX(t) = X; tt(t) = t; yy(t) = y;
    
end



%% Display Results

figure(); ...
    plot(tt, xx, 'k', 'Marker', 'o'); hold on;
    plot(tt, xx);
    plot(tt, xx_predicted);
    plot(tt, yy, 'b', 'Marker', 's');
    legend('Ground Truth', 'Corrected State', 'Predicted State', 'Measurements', 'Location', 'SouthEast');
    xlabel('Time (seconds)'); ylabel('Position'); title('1D Kalman Filter');
    grid on; shg;
    
    
%% Clean-up

fprintf(1, '\n*** Processing Complete ***\n');



%% Reference(s)


