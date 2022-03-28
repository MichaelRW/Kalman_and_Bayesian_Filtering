

%% Synopsis:

% Example using a Kalman filter to
% estimate the resting membrane potential of a neuron
% (Algorithm isn't optimized)



%% Environment

close all; clear; clc;

set( 0, 'DefaultFigureWindowStyle', 'docked' );



%% Parameters

N=1000;         % number of measurements
Vm=-70;         % the membrane potential (Vm) value
sw=1^2;         % variance process noise
sv=2^2;         % measurement noise
% Initial Estimates
x_apriori(1)=Vm;        % first estimate of Vm
                        % This is the a priori measurement
s_apriori(1)=0;         % a priori estimate of the error; here variance =0 
                        % i.e., we assume that first estimate is perfect

% create a simulated measurement vector
% -------------------------------------
for i=1:N;
    tru(i)=Vm+randn*sw;
    z(i)=tru(i)+randn*sv;
end;

% Perform the Kalman filter
% -------------------------
% Note: Malab indices start at 1, and not
%       at 0 as in most textbook derivations

% Equations are numbered as in CH 19
for i=1:N;
    % step 1: Blending of Prediction and Measurement 
    % Compute blending factor K(i)
    K(i)=s_apriori(i)/(s_apriori(i)+sv);                      % Eq (19.19)
    % Update Estimate using the measurement
    % this is the a posteriori estimate
    x_aposteriori(i)=x_apriori(i)+K(i)*(z(i)-x_apriori(i));   % Eq (19.15)
    % Update the error estimate [Note that there
    % are several variants for this procedure;
    % here we use the simplest expression
    s_aposteriori(i)=s_apriori(i)*(1-K(i));                   % Eq (19.20)
    % step 2: Project Estimates for the next Step
    x_apriori(i+1)=x_aposteriori(i);                          % Eq (19.21)
    s_apriori(i+1)=s_aposteriori(i)+sw;                       % Eq (19.22)
end;



%% Plot Results

figure; ...
    plot( z, 'k' );  hold on;
    plot( tru, '-.r' );
    plot( x_aposteriori, '--b' );  grid on;
    axis( [ 3e2 3.5e2 -85 -55 ] );
    legend( 'True Values' , 'Measurements' ,'Kalman Filter Output' );
    title ( 'Application - Scalar Kalman Filter' );
    xlabel ( 'Time (AU) ' );  ylabel ( 'Trans-membrane Potential (mV)' );
    %
    shg

fprintf( 1, '\nVariance - Truth:  %4.2f\tMeasurement:  %4.2f\tx posteriori:  %4.2f\n\n', var(tru), var(z), var(x_aposteriori) );



%% Clean-up

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)


