

%% Acknowledgements

% Written by:  StudentDave
%
% For licensing and usage questions e-mail:  scienceguy5000 at gmail. com

% Recursive Bayesian estimation example adapted from Michael A. Goodrich - September 3, 2004 - CS 470



%% Environment

close all; clear; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');

addpath(genpath('.'), '-begin');

h1=figure('Name', 'Quail and Apparent Squawk Locations');
h2=figure('Name', 'Posterior PMF');
h3=figure('Name', 'Position Updates');
    autoArrangeFigures();



%% Simulation Parameters

numberOfHeardSquawks=50;
quailLocation=[ 5; 5 ];

% The hidden quail squawks "numberOfHeardSquawks".  For each squawk, the ninja hears the quail
% at a random location.  The likelihood of this apparent location and the physical location of the
% quail will be estimate as a 2-dimensional Gaussian distribution centered at the quail's location.
% The standard deviation of the 2-dimensional Gaussion distribution is 2 in both dimensions.

n=2*randn(2, numberOfHeardSquawks);  % 2-dimensional Gaussian squawk distribution.
apparentQuailLocations=zeros(2, numberOfHeardSquawks);  % Apparent location of the quail (i.e. the source of the squawk).



%% Display Squawk Locations

figure(h1); ...
    plot( quailLocation(1), quailLocation(2), 'r.', 'MarkerSize', 40, 'LineWidth', 3 );  % Location of quail.
    hold on; for i=1:numberOfHeardSquawks; apparentQuailLocations(:, i)=quailLocation + n(:, i); plot(apparentQuailLocations(1, i), apparentQuailLocations(2, i), 'k.', 'markersize', 10); end;
    axis( [-5, 15, -5, 15] ); daspect([1 1 1]); grid on; legend('Location of Quail', 'Apparent Squawk Locations');
    
    
    
%% Set Quail Locations, Prior and Posterior PMFs.

xLocations=4:0.05:6; yLocations=xLocations;  % Define the locations, or states, where the quail can potentially be located.

% Set the prior knowledge about the quail's location.
PRIOR_KNOWLEDGE=3;

switch PRIOR_KNOWLEDGE
    case 1
        % No bias, no prior knowledge, and uniform distribution.
        L=length(xLocations); priorStateKnowledge=ones(L); posteriorStateKnowledge=ones(L);
        
    case 2
        % Bias the prior towards the quail's location.
        %
        % Generated Gaussian distribution centered at given point.
        %
        centered_prior=zeros(size(xLocations, 2), size(yLocations, 2));
        
        for xIndex=1:1:numel(xLocations)            
            for yIndex=1:1:numel(yLocations)
                centered_prior(xIndex, yIndex)=exp( -1/2 .* ([xLocations(xIndex); yLocations(yIndex)] - [5; 5]).' * ( [0.05 0; 0 0.05] \ ([xLocations(xIndex); yLocations(yIndex)] - [5; 5]) ) ) / ...
                    sqrt( (2*pi)^2 * det([0.05 0; 0 0.05]) );
            end            
        end
        
        priorStateKnowledge=centered_prior; posteriorStateKnowledge=centered_prior;
            
    case 3
        % Bias the prior away from the quail's location.
        %
        % Generated Gaussian distribution centered at given point.
        %
        offCentered_prior=zeros(size(xLocations, 2), size(yLocations, 2));
        
        for xIndex=1:1:numel(xLocations)            
            for yIndex=1:1:numel(yLocations)
                offCentered_prior(xIndex, yIndex)=exp( -1/2 .* ([xLocations(xIndex); yLocations(yIndex)] - [4; 4]).' * ( [0.05 0; 0 0.05] \ ([xLocations(xIndex); yLocations(yIndex)] - [4; 4]) ) ) / ...
                    sqrt( (2*pi)^2 * det([0.05 0; 0 0.05]) );
            end            
        end
        
        priorStateKnowledge=offCentered_prior; posteriorStateKnowledge=offCentered_prior;
                    
    otherwise
        error( '*** Invalid prior knowledge CASE index. ***' );
end;


% Convert the state knowledge into a probability mass function (PMF; discrete values).
priorStateKnowledge=priorStateKnowledge / sum(priorStateKnowledge(:));
posteriorStateKnowledge=posteriorStateKnowledge / sum(posteriorStateKnowledge(:));

% figure; mesh(priorStateKnowledge); shg;
% figure; mesh(posteriorStateKnowledge); shg; return;


% Display the posterior knowledge matrix (see SWITCH argument above).
figure(h2); clf; ...
    mesh(posteriorStateKnowledge); colorbar; ...
    xticklabels(4:0.5:6); yticklabels(4:0.5:6); 
    axis([0 40 0 40 0 0.05]); ...
    xlabel('X'); ylabel('Y'); zlabel('Probability');



%% Iterative Bayes Computations

% Find indices where the posteriorStateKnowledge has its maximum value.
[maximumX_posteriorPMF_index, maximumY_posteriorPMF_index]=find( posteriorStateKnowledge==max(posteriorStateKnowledge(:)) );
%
% Note:  With the three cases above, the maximum x and y PMF value indices will be,
%   Case 1 - With a uniform posterior PMF (i.e. from the prior PMF) there will not be an optimal state value.
%   Case 2 - A Gaussian posterior PMF with a maximum value near the quail's position.
%   Case 3 - A Gaussion posterior PMF with a maximum value offset from the quail's position.

stateEstimate=[ xLocations(maximumX_posteriorPMF_index); yLocations(maximumY_posteriorPMF_index) ];

figure(h3); ...
    subplot(2,1,1); ...
        plot(1, stateEstimate(1)); hold on;
        line([1, numberOfHeardSquawks], [quailLocation(1), quailLocation(1)] );
    subplot(2,1,2); ...
        plot(1, stateEstimate(2)); hold on;
        line([1, numberOfHeardSquawks], [quailLocation(2), quailLocation(2)] );
        
twoDimensionCovarianceMatrix=[ 4, 0; 0, 4 ];

for dataPointIndex=2:1:length(apparentQuailLocations)
    
    priorStateKnowledge=posteriorStateKnowledge;
        temporaryPosteriorStateKnowledge=0*priorStateKnowledge;  % FIXME
        
    % Compute the likelihood that the apparent location of the quail's squawk was made
    % at each possible location in the bush.
    %
    for xPositionIndex=1:1:length(priorStateKnowledge)
        
        for yPositionIndex=1:1:length(priorStateKnowledge)
            
            workingPosition=[ xLocations(xPositionIndex); yLocations(yPositionIndex) ];
            
            gaussianDenominator=sqrt( (2*pi)^2 * det(twoDimensionCovarianceMatrix) );
            
            temporaryPosteriorStateKnowledge(xPositionIndex, yPositionIndex)=exp( -1/2 * ...
                ( apparentQuailLocations(:, dataPointIndex) - workingPosition ).' * ...
                ( twoDimensionCovarianceMatrix \ (apparentQuailLocations(:, dataPointIndex) - workingPosition) ) ) / gaussianDenominator;
            
            temporaryPosteriorStateKnowledge(xPositionIndex, yPositionIndex)=temporaryPosteriorStateKnowledge(xPositionIndex, yPositionIndex) * priorStateKnowledge(xPositionIndex, yPositionIndex);
            
        end;
        
    end;
    
    posteriorStateKnowledge=temporaryPosteriorStateKnowledge / sum(temporaryPosteriorStateKnowledge(:));
    
    
    figure(h2); ...
        mesh(posteriorStateKnowledge); colorbar; ...
        xticklabels(4:0.5:6); yticklabels(4:0.5:6); 
        axis([0 40 0 40 0 0.05]); ...
        xlabel('Y'); ylabel('X'); zlabel('Probability');
    
    [maximumX_posteriorPMF_index, maximumY_posteriorPMF_index]=find( posteriorStateKnowledge==max(max(posteriorStateKnowledge)) );
    stateEstimate=[ xLocations(maximumX_posteriorPMF_index); yLocations(maximumY_posteriorPMF_index) ];
    
    figure(h3); ...    
        subplot(2,1,1); ...
            plot(dataPointIndex, stateEstimate(1), 'k.'); axis([ 0 numberOfHeardSquawks 4 6 ]); ylabel('X Estimate'); grid on;
        subplot(2,1,2); ...
            plot(dataPointIndex, stateEstimate(2), 'k.'); axis([ 0 numberOfHeardSquawks 4 6 ]); xlabel('Data Index'); ylabel('Y Estimate'); grid on;

end;

fprintf(1, '\n\nFinal location estimate:  X: %3.1f, Y: %3.1f\n\n', stateEstimate(1), stateEstimate(2));



%% Reference(s)

% https://en.wikipedia.org/wiki/Multivariate_normal_distribution


