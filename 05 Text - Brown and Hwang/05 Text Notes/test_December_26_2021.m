

% https://www.youtube.com/watch?v=4jRBRDbJemM&t=275s&ab_channel=StatQuestwithJoshStarmer


%% Environment

close all force; clc; clear; restoredefaultpath;
set( 0, 'DefaultFigureWindowStyle', 'normal' );



%% White Noise

% AVERAGING_SIZE = 10;
%     b = ones(1, AVERAGING_SIZE) ./ AVERAGING_SIZE;  a = 1;
    b = 1; a = [ 1 -0.5 ];
    %
    figure; freqz( b, a, 4096 );

rng( 0 );
    x = randn( 1, 2^15 );
        xF = filter( b, a, x );
        
corrAt0 = corr2( xF(1:end).', xF(1:end).' )
corrAt1 = corr2( xF(1:end-1).', xF(2:end).' )
corrAt2 = corr2( xF(1:end-2).', xF(3:end).' )
    
LAGS = 20;        
    [ r, lags ] = xcorr( x, LAGS, 'coeff' );
    [ rF, lagsF ] = xcorr( xF, LAGS, 'coeff' );
%
figure; ...
    plot( lags, r, 'Marker', '.' ); hold on;
    stem( lagsF, rF, 'r', 'Marker', 'o' ); grid on; shg;
    axis( [ -LAGS-5 LAGS+5, -0.15, 1.5] );
    legend( 'White Noise', 'Coloured Noise' );
    
    
    
%% Clean-up

autoArrangeFigures(2, 2, 1);

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );


