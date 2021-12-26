

% https://www.youtube.com/watch?v=4jRBRDbJemM&t=275s&ab_channel=StatQuestwithJoshStarmer


%% Environment

close all force; clc; clear; restoredefaultpath;
set( 0, 'DefaultFigureWindowStyle', 'normal' );



%% Deterministic Random Process

rng( 0 );  VARIANCE = 2;
    A = sqrt( VARIANCE ) .* randn( 50, 1 );
    
netTime_seconds = 5;  sampleRate_Hz = 1e3;  samplePeriod_seconds = 1/sampleRate_Hz;
    timeIndices = samplePeriod_seconds:samplePeriod_seconds:(netTime_seconds - samplePeriod_seconds);
    
x = A * sin( 2 * pi * 20 .* timeIndices );
%
% figure; plot( x.' ); grid on; shg;

% Time Autocorrelation
corrAt7 = corr2( x(1, 1:end-7).', x(1, 8:end).' )

LAGS = 50;        
    [ r, lags ] = xcorr( x(1, :), LAGS, 'coeff' );
%
figure; plot( lags, r, 'Marker', '.' ); grid on; shg;


% Ensemble Autocorrelation
VARIANCE .* sin( 2 * pi * 20 * samplePeriod_seconds ) .* sin( 2 * pi * 7 * samplePeriod_seconds )
% ensembleSet = A .* sin( 2 * pi * 20 .* timeIndices );

tempA = VARIANCE .* sin( 2 * pi * 20 * timeIndices ) .* sin( 2 * pi * 20 * (timeIndices + 7 * samplePeriod_seconds) );
    figure; plot( tempA ); grid on; shg;
    
    
    
%% Clean-up

autoArrangeFigures(2, 2, 1);

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );


