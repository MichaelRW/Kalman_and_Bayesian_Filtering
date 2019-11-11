

%% Environment

close all; clear; clc;
set(0, 'DefaultFigureWindowStyle', 'docked');
rng(0);



%% Entities

NUMBER_OF_MEASUREMENTS = 1e4;

nominal_resistance = 100; resistors_rms_error = 1;
    resistors = nominal_resistance + randn(NUMBER_OF_MEASUREMENTS, 1)*resistors_rms_error;

measurements_rms_error = 3;
    measurements = nominal_resistance + randn(NUMBER_OF_MEASUREMENTS, 1)*measurements_rms_error;
        figure; stem(measurements); title('Measurements'); grid on; shg;


%% Least Squares Estimation

leastSquares_RMS_error = measurements_rms_error ./ (1:1e-1:NUMBER_OF_MEASUREMENTS);
    h = figure; ...
        plot(1:1e-1:NUMBER_OF_MEASUREMENTS, leastSquares_RMS_error);
            ylabel('RMS Error'); xlabel('Observation'); title('Estimate RMS Error');
            xlim([0 NUMBER_OF_MEASUREMENTS]); grid on; shg;
            
            
            
%% Maximum Likelihood Estimation

OHM_METER_SCALE_FACTOR = 1;

b = zeros(NUMBER_OF_MEASUREMENTS, 1);
    b(1) = (OHM_METER_SCALE_FACTOR * resistors_rms_error^2 ) / (OHM_METER_SCALE_FACTOR^2 * resistors_rms_error^2 + measurements_rms_error^2);
    
estimates = nan(NUMBER_OF_MEASUREMENTS, 1); errorVariance = nan(NUMBER_OF_MEASUREMENTS, 1);

estimates(1) = nominal_resistance + b(1)*(measurements(1) - nominal_resistance);
errorVariance(1) = (1 - b(1)*OHM_METER_SCALE_FACTOR)*resistors_rms_error^2;

for estimateIndex = 2:1:NUMBER_OF_MEASUREMENTS
        
    b(estimateIndex) = (OHM_METER_SCALE_FACTOR * errorVariance(estimateIndex - 1)) / (OHM_METER_SCALE_FACTOR^2 * errorVariance(estimateIndex - 1) + measurements_rms_error^2);
    
    estimates(estimateIndex) = estimates(estimateIndex - 1) + b(estimateIndex)*(measurements(estimateIndex) - OHM_METER_SCALE_FACTOR*estimates(estimateIndex - 1));
    
    errorVariance(estimateIndex) = (1 - b(estimateIndex)*OHM_METER_SCALE_FACTOR)*errorVariance(estimateIndex - 1);
    
end
    
figure; ...
    subplot(1,3,1); plot(b); grid on; title('Weighting (b)');
    subplot(1,3,2); plot(estimates, 'b'); grid on; title('Estimate');
    subplot(1,3,3); plot(errorVariance, 'r'); grid on; title('Error Variance');
    
    
    