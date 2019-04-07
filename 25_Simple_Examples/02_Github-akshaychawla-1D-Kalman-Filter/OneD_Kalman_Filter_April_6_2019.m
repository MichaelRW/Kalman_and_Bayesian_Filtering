

%% Environment

close all; clear; clc;
set(0, 'DefaultFigureWindowStyle', 'docked');



%% System Variables

X=5; % Initial position state estimate.
timeStep_seconds=1;  % State update time interval.
velocity_m_per_s=5; %speed/control = 2m/s

rng(0);  % Fix pseudo-random number generator to produce identical values.


%% Kalman Filter Variables

% Initial State Values (1D System)
x=0;  % Initial position of vehicle.
P=1;  % Variance estimate of initial position.

% Time Update Parameters
A=1;  % State transition matrix.
B=1;  % Control input matrix.
Q=1;  % Process noise variance.

% Measurement Update Parameters
H=1;  % Measurement matrix.
R=4;  % Measurement noise variance.



%% Initialize Processing Variables

simulationSize=1e2;
    XX=zeros(1, simulationSize);
    tt=zeros(1, simulationSize);
    xx_predicted=zeros(1, simulationSize);
    xx=zeros(1, simulationSize);
    PP=zeros(1, simulationSize);
    yy=zeros(1, simulationSize);
    kk=zeros(1, simulationSize);



%% Simulating the System

% h1=figure('Name', 'Time Update Gaussian');

for t=1:1:simulationSize
    
    n=sqrt(Q) * randn();  % State prediction noise (zero mean Gaussian distribution).
        X=X + velocity_m_per_s*timeStep_seconds + n;  % Predict next state.
        
    v=sqrt(R) * randn();  % State measurement noise (zero mean Gaussian distribution).
        y=H*X + v;  % Measurement of state.
    
        
   % Time Update (Prediction) -  Prior 
    x_predicted=A*x + B*velocity_m_per_s;  % Predicting the new state (mean).
        P=A*P*A.' + Q;  % Update the uncertainty of the new state.
    %
%     figure(h1); hold on; plot([0:0.1:200], normpdf([0:0.1:200], x_predicted, sqrt(P))); grid on; shg;
    
    
    % Measurement Update (Correction) -  Posterior
    K=(H*P*H.'+R)\P*H.';  % Compute Kalman gain.    
        x=x_predicted+K*(y-H*x_predicted);  % Update the estimate via y.    
        P=(1-K*H)*P;
	%
%     figure(h1); hold on; plot([0:0.1:200], normpdf([0:0.1:200], x, sqrt(P)), 'r'); shg; pause(1);    
   

	% Store data for plotting prior state values.
    xx_predicted(t)=x_predicted;
    xx(t)=x; PP(t)=P; XX(t)=X; tt(t)=t; yy(t)=y; kk(t)=K;
    
end



%% Display Results

figure(); ...
    plot(tt, XX, 'k', 'Marker', 'o'); hold on;
    plot(tt, xx, 'r', 'Marker', '^');
    plot(tt, xx_predicted);
    plot(tt, yy, 'b', 'Marker', 's');
    legend('Ground Truth', 'Corrected State (a Posterior)', 'Predicted State (a Prior)', 'Measurements', 'Location', 'SouthEast');
    xlabel('Time (seconds)'); ylabel('Position'); title('1D Kalman Filter');
    grid on; shg;
    
figure(); plot(tt, PP); grid on; shg;

figure(); plot(tt, kk); grid on; shg;

figure(); plot(abs(XX-xx)); hold on; plot(abs(XX-xx_predicted), 'r'); grid on; shg; legend('a Posterior', 'a Prior');

    
    
%% Clean-up

fprintf(1, '\n*** Processing Complete ***\n');



%% Reference(s)

% https://www.mathworks.com/examples/statistics/mw/stats-ex15608473-compute-and-plot-the-normal-distribution-pdf


