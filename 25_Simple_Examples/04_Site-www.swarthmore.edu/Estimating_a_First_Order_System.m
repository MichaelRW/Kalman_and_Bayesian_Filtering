

%% Synopsis

% 



%% Environment

close all; clear; clc;
set(0, 'DefaultFigureWindowStyle', 'docked');



%% Simulation Parameters

stateIndices=25;  % Length of the simulation.

stateTransitionMatrix=0.85;  % 1 for a constant; |a|<1 for a first-order system.
    Q=0.0;  % Process noise variance.

measurementMatrix=3;
    R=1;  % Measurement noise variance.

% Preallocated memory for result vectors.
x=zeros(1, stateIndices);
z=zeros(1, stateIndices);
xapriori=zeros(1, stateIndices);
xaposteriori=zeros(1, stateIndices);
residual=zeros(1, stateIndices);
papriori=ones(1, stateIndices);
paposteriori=ones(1, stateIndices);
k=zeros(1, stateIndices);

% Preallocate the process and measurement noise.
rng(0);
    w1=randn(1, stateIndices); w=w1.*sqrt(Q);
    v1=randn(1, stateIndices); v=v1.*sqrt(R);
    
    
    
%% Processing - Intialization

% Initial condition on the position state, x.
x_0=1.0;

% Initial guesses for state and a posteriori covariance.
xaposteriori_0=1.5;
paposteriori_0=1;

% Calculate the first estimates for all values based upon the initial guess
% of the state and the a posteriori covariance.  The rest of the steps will
% be calculated in a loop.

% Calculate the state and the output.
x(1)=stateTransitionMatrix*x_0 + w(1);
    z(1)=measurementMatrix*x(1) + v(1);
    
% Predictor equations.
xapriori(1)=stateTransitionMatrix*xaposteriori_0;
    residual(1)=z(1) - measurementMatrix*xapriori(1);
papriori(1)=stateTransitionMatrix*stateTransitionMatrix*paposteriori_0 + Q;

% Corrector equations.
k(1)=measurementMatrix*papriori(1) / (measurementMatrix*measurementMatrix*papriori(1) + R);
    paposteriori(1)=papriori(1)*(1 - measurementMatrix*k(1));
    xaposteriori(1)=xapriori(1) + k(1)*residual(1);



%% Processing - Iterations/Recursions

for stateIndex=2:stateIndices
    
    % Calculate the state and the output.
    x(stateIndex)=stateTransitionMatrix*x(stateIndex-1) + w(stateIndex);
        z(stateIndex)=measurementMatrix*x(stateIndex) + v(stateIndex);
    
    % Predictor equations.
    xapriori(stateIndex)=stateTransitionMatrix*xaposteriori(stateIndex-1);
        residual(stateIndex)=z(stateIndex) - measurementMatrix*xapriori(stateIndex);
    papriori(stateIndex)=stateTransitionMatrix*stateTransitionMatrix*paposteriori(stateIndex-1)+Q;
    
    % Corrector equations.
    k(stateIndex)=measurementMatrix*papriori(stateIndex) / (measurementMatrix*measurementMatrix*papriori(stateIndex)+R);
        paposteriori(stateIndex)=papriori(stateIndex)*(1 - measurementMatrix*k(stateIndex));
        xaposteriori(stateIndex)=xapriori(stateIndex) + k(stateIndex)*residual(stateIndex);
    
end



%% Display Results

stateIndex=1:stateIndices;

figure(); ...
    
    subplot(2,2,1); ...
        h1=plot(stateIndex, xapriori, 'b', 'Marker', 'v'); hold on;
        h2=plot(stateIndex, xaposteriori, 'g', 'Marker', '^');
        h3=plot(stateIndex, x, 'r', 'Marker', '.'); hold off;
        legend([h1(1) h2(1) h3(1)], 'a priori', 'a posteriori', 'exact', 'Location', 'NorthEast');
        title('State with a priori and a posteriori Estimate'); ylabel('State, x');
        xlim=[0 length(stateIndex)+1]; set(gca,'XLim', xlim); grid on;
        
    subplot(2,2,2); ...
        h1=plot(stateIndex, papriori, 'b', 'Marker', 'v'); hold on;
        h2=plot(stateIndex, paposteriori, 'g', 'Marker', '^'); hold off
        legend([h1(1) h2(1)], 'a priori', 'a posteriori');
        title('Calculated a priori and a posteriori Variance'); ylabel('Variance');
        set(gca,'XLim', xlim); grid on;
        
	subplot(2,2,3); ...
        h1=stem(stateIndex, x-xapriori, 'b', 'Marker', 'v'); hold on;
        h2=stem(stateIndex, x-xaposteriori, 'g', 'Marker', '^'); hold off;
        legend([h1(1) h2(1)], 'a priori', 'a posteriori', 'Location', 'SouthEast');
        title('Actual a priori and a posteriori error'); ylabel('Errors');
        set(gca,'XLim',xlim); grid on;
        
    subplot(2, 2, 4); ...
        h1=stem(stateIndex, k, 'b');
        legend(h1(1), 'Kalman Gain');
        title('Kalman Gain (K)'); ylabel('Kalman Gain (K)');
        set(gca,'XLim',xlim); grid on;
        
        shg;
        
        
        
%% Clean-up

fprintf(1, '\n\n*** Processing Complete ***\n\n');



%% Reference(s)

% http://www.swarthmore.edu/NatSci/echeeve1/Ref/Kalman/scalarkalman.m


