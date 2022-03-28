% pr19_1
%
% Example using a Kalman filter to
% estimate the resting membrane potential of a neuron
% (Algorithm isn't optimized)


clear
close all

% Parameters
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

% plot the results
figure;hold;
plot(z,'r.');
plot(tru,'g');
plot(x_aposteriori)
axis([1 N -100 10]);
xlabel ('Time (AU)')
ylabel ('Membrane Potential (mV)')
title ('Kalman Filter Application: true values (green); measurements (red .); Kalman Filter Output (blue)')
