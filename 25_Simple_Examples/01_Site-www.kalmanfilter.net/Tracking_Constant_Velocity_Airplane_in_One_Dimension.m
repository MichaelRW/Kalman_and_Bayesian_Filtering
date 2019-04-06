

%% Environment

close all; clear; clc;
set(0, 'DefaultFigureWindowStyle', 'docked');



%% Example from Website

alpha=0.2; beta=0.1; deltaT_seconds=5;
% alpha=0.8; beta=0.8; deltaT_seconds=5;

initialPosition_meters=30e3; initialVelocity_meterPerSecond=40;

predictedStateEstimate=zeros(2, 12);  % Row 1 - Position, Row 2 - Velocity
extrapolatedStateEstimate=zeros(2, 12);
test=zeros(2, 12);

distanceMeasurements_meters=[ 30110, 30265 30740 30750 31135 31015 31180 31610 31960 31865 ];

predictedStateEstimate(1, 1)=initialPosition_meters + deltaT_seconds*initialVelocity_meterPerSecond;
predictedStateEstimate(2, 1)=initialVelocity_meterPerSecond;


for stateUpdateIndex=2:1:10
    
    % Position
    predictedStateEstimate(1, stateUpdateIndex)=predictedStateEstimate(1, stateUpdateIndex-1) + ...
        alpha*(distanceMeasurements_meters(stateUpdateIndex-1) - predictedStateEstimate(1, stateUpdateIndex-1));
    %
    % Velocity
    predictedStateEstimate(2, stateUpdateIndex)=predictedStateEstimate(2, stateUpdateIndex-1) + ...
        beta*( (distanceMeasurements_meters(stateUpdateIndex-1) - predictedStateEstimate(1, stateUpdateIndex-1) ) / deltaT_seconds );
    
    
    % Position
    extrapolatedStateEstimate(1, stateUpdateIndex)=predictedStateEstimate(1, stateUpdateIndex) + ...
        deltaT_seconds*predictedStateEstimate(2, stateUpdateIndex);
    %
    % Velocity
    extrapolatedStateEstimate(2, stateUpdateIndex)=predictedStateEstimate(2, stateUpdateIndex);
    
    
    test(1, stateUpdateIndex-1)=predictedStateEstimate(1, stateUpdateIndex-1) + ...
        alpha*(distanceMeasurements_meters(stateUpdateIndex-1) - predictedStateEstimate(1, stateUpdateIndex-1));
    %
    test(2, stateUpdateIndex-1)=predictedStateEstimate(2, stateUpdateIndex-1) + ...
        beta*((distanceMeasurements_meters(stateUpdateIndex-1) - predictedStateEstimate(1, stateUpdateIndex-1))/deltaT_seconds);
    
    
    predictedStateEstimate(1, stateUpdateIndex)=extrapolatedStateEstimate(1, stateUpdateIndex);
    predictedStateEstimate(2, stateUpdateIndex)=extrapolatedStateEstimate(2, stateUpdateIndex);
    
%     keyboard;
    
end


figure(); ...   
    plot((30200+40*(0:1:9)*deltaT_seconds), 'g'); hold on;
    plot(distanceMeasurements_meters, 'b', 'Marker', 's');
    plot(predictedStateEstimate(1, 1:10), 'k', 'Marker', '^');
    plot(1:9, test(1, 1:9), 'r', 'Marker', 'o');
    xlabel('Measurement Index'); ylabel('Range (kilometers)'); title('Range Estimate - 1D Kalman Processing');
    legend('True Value', 'Measurements', 'Prediction', 'Estimates', 'Location', 'SouthEast');
    grid on; shg;



%% Clean-up

fprintf(1, '\n\n*** Processing Complete ***\n\n');



%% References


