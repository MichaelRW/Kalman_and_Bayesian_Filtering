 

%% Environment

close all; clear; clc;

% addpath(genpath(fullfile('..', '..', '95 Reference')), '-begin');

set(0, 'DefaultFigureWindowStyle', 'docked');



%% Example from Website

initialWeightEstimate_grams=1.03e3;

weightMeasurements_grams=[1.03 0.989 1.017 1.009 1.013 0.979 1.008 1.042 1.012 1.011].*1e3;

gain=1./(1:1:numel(weightMeasurements_grams));

priorStateEstimate=zeros(1, numel(weightMeasurements_grams)); priorStateEstimate(1)=initialWeightEstimate_grams;
    
updatedStateEstimate=zeros(1, numel(weightMeasurements_grams));

trueStateValue_grams=1.01e3*ones(1, numel(weightMeasurements_grams));


for stateUpdateIndex=2:1:10    
    updatedStateEstimate(stateUpdateIndex) = priorStateEstimate(stateUpdateIndex-1) + ...
        gain(stateUpdateIndex)*( weightMeasurements_grams(stateUpdateIndex) - priorStateEstimate(stateUpdateIndex-1) );
    %
    priorStateEstimate(stateUpdateIndex) = updatedStateEstimate(stateUpdateIndex);
end


figure(); ...
    plot(trueStateValue_grams, 'g'); hold on;
    plot(weightMeasurements_grams, 'b', 'Marker', 's');
    plot(priorStateEstimate, 'r', 'Marker', 'o');
    xlabel('Measurement Index'); ylabel('Weight (grams)'); title('Weight Estimate - 1D Kalman Processing');
    legend('True Value', 'Measurements', 'Estimates', 'Location', 'SouthEast');
    grid on; shg;



%% Clean-up

fprintf(1, '\n\n*** Processing Complete ***\n\n');



%% References


