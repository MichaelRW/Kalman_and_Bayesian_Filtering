

%% Environment

close all; clear; clc;
set(0, 'DefaultFigureWindowStyle', 'docked');



%% System Variables

X=0; % Initial position state estimate.
timeStep_ms=1;  % State update time interval.

rng(0);  % Fix pseudo-random number generator to produce identical values.


%% Kalman Filter Variables

x=0;  % Initial position of vehicle.
P=1;  % Variance estimate of initial position.

% Time Update Parameters
A=1;  % State transition matrix.
B=0;  % Control input matrix.
Q=0.01;  % Process noise variance.

% Measurement Update Parameters
H=1;  % Measurement matrix.
R=0.1;  % Measurement noise variance.

% measurements=[ 0.39 0.50 0.48 0.29 0.25 0.32 0.34 0.48 0.41 0.45];
measurements=sqrt(0.5)*randn(1, 1e2)+0.5;



%% Initialize Processing Variables

simulationSize=numel(measurements);
    XX=zeros(1, simulationSize);
    tt=zeros(1, simulationSize);
    xx_predicted=zeros(1, simulationSize);
    xx=zeros(1, simulationSize);
    PP=zeros(1, simulationSize);
    yy=zeros(1, simulationSize);
    kk=zeros(1, simulationSize);



%% Simulating the System

% h1=figure('Name', 'Time Update Gaussian');

for t=1:1:numel(measurements);
    
   % Time Update (Prediction) -  Prior 
    x_predicted=A*x;  % Predicting the new state (mean).
        P=A*P*A.' + Q;  % Update the uncertainty of the new state.
    %
%     figure(h1); hold on; plot([0:0.1:200], normpdf([0:0.1:200], x_predicted, sqrt(P))); grid on; shg;
    
    
    % Measurement Update (Correction) -  Posterior
    K=(H*P*H.'+R)\P*H.';  % Compute Kalman gain.    
        x=x_predicted+K*(measurements(t)-H*x_predicted);  % Update the estimate via y.    
        P=(1-K*H)*P;
	%
%     figure(h1); hold on; plot([0:0.1:200], normpdf([0:0.1:200], x, sqrt(P)), 'r'); shg; pause(1);    
   

	% Store data for plotting prior state values.
    xx_predicted(t)=x_predicted;
    xx(t)=x; PP(t)=P; XX(t)=X; tt(t)=t; kk(t)=K;
    
end



%% Display Results

figure(); ...
    plot(tt, xx, 'k', 'Marker', 'o'); hold on;
    plot(tt, xx_predicted, 'r', 'Marker', '^');
    plot(tt, measurements, 'b', 'Marker', 's');
    plot(tt, 0.5*ones(1, numel(tt)), 'm', 'Marker', 'v');
    legend('Predicted (a Prior)', 'Measurement Update (a Posterior)', 'Measurements', 'Mean Voltage', 'Location', 'SouthEast');
    xlabel('Time (seconds)'); ylabel('Volts'); title('1D Kalman Filter');
    grid on; shg;
    
figure(); plot(tt, PP); grid on; shg;

figure(); plot(tt, kk); grid on; shg;

    
    
%% Clean-up

fprintf(1, '\n*** Processing Complete ***\n');



%% Reference(s)

% http://bilgin.esme.org/BitsAndBytes/KalmanFilterforDummies#

% https://www.mathworks.com/examples/statistics/mw/stats-ex15608473-compute-and-plot-the-normal-distribution-pdf


