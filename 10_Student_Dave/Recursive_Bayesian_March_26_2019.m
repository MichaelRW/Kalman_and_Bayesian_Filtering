

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
h3=figure(3);
    autoArrangeFigures();



%% Simulation Parameters

numberOfHeardSquawks=50;
quailLocation=[ 5; 5 ];

% The hidden quail squawks "numberOfHeardSquawks".  For each squawk, the ninja hears the quail
% at a random location.  The likelihood of this apparent location and the physical location of the
% quail will be estimate as a 2-dimensional Gaussian distribution centered at the quail's location.
% The standard deviation of the 2-dimensional Gaussion distribution is 2 in both dimensions.

n=2*randn(2, numberOfHeardSquawks);  % 2-dimensional Gaussian squawk distribution.
x=zeros(2, numberOfHeardSquawks);  % Apparent location of the quail (i.e. the source of the squawk).



%% Display Squawk Locations

figure(h1); ...
    plot( quailLocation(1), quailLocation(2), 'r.', 'MarkerSize', 40, 'LineWidth', 3 );  % Location of quail.
    hold on; for i=1:numberOfHeardSquawks; x(:, i)=quailLocation + n(:, i); plot(x(1, i), x(2, i), 'k.', 'markersize', 10); end;
    axis( [-5, 15, -5, 15] ); daspect([1 1 1]); grid on; legend('Location of Quail', 'Apparent Squawk Locations');
    
    
    
%% Bayesian Analysis

xLocations=4:0.05:6; yLocations=xLocations;  % Define the locations, or states, where the quail can potentially be located.

% Set the prior knowledge about the quail's location.
PRIOR_KNOWLEDGE=1;

switch PRIOR_KNOWLEDGE
    case 1
        % No bias, no prior knowledge, and uniform distribution.
        L=length(xLocations); priorStateKnowledge=ones(L, L); posteriorStateKnowledge=ones(L, L);
    case 2
        % Bias the prior towards the quail's location.
        centered_prior=[4, 4];
            priorStateKnowledge=centered_prior; posteriorStateKnowledge=centered_prior;
    case 3
        % Bias the prior away from the quail's location.
        off_centered_prior=[0, 0];
            priorStateKnowledge=off_centered_prior; posteriorStateKnowledge=off_centered_prior;
    otherwise
        error( '*** Invalid prior knowledge CASE index. ***' );
end;


% Convert the state knowledge into a probability mass function (PMF; discrete values).
priorStateKnowledge=priorStateKnowledge / sum(priorStateKnowledge(:));
posteriorStateKnowledge=posteriorStateKnowledge / sum(posteriorStateKnowledge(:));

% Display the posterior knowledge matrix;
figure(h2); clf; mesh(posteriorStateKnowledge); xticklabels(4:0.5:6); colorbar; axis([0 40 0 40 0 0.015]); daspect([1 1 1]);

return;




%%iterative bayes

[a,b]=find(posteriorStateKnowledge==max(max(posteriorStateKnowledge)));  % Pull out the indices at which posteriorStateKnowledge achieves its max to start.
sest=[xLocations(a);yLocations(b)];  % The best estimate of the true state to start.
figure(h1);
clf
figure(h2);
clf
subplot(211); plot(1,sest(1)); hold on;
line([1,numberOfHeardSquawks],[quailLocation(1),quailLocation(1)]); % Draw a line at the location of the x dimension.
subplot(212); plot(1,sest(2)); hold on;
line([1,numberOfHeardSquawks],[quailLocation(2),quailLocation(2)]); % Draw a line at the location of the y dimension.

K=[4,0;0,4]; % covariance matrix for making a 2-D gaussian
for (n=2:length(x));
    priorStateKnowledge=posteriorStateKnowledge; %store the posterior to the prior.
    m=0*priorStateKnowledge;   
    %likelihood
    % look at each location, assume that the given location is
    % is there quail is, and get the likelihood of the data x(:,n) assuming
    % 2-d gaussian noise
    for (i=1:length(priorStateKnowledge))
       for (j=1:length(priorStateKnowledge))
           me=[xLocations(i); yLocations(j)];
           m(i,j) = 1/sqrt((2*pi)^2*det(K)) * exp(-(x(:,n)-me)'*inv(K)*(x(:,n)-me)/2); %Compute likelihood           
           m(i,j) = m(i,j) * priorStateKnowledge(i,j); % Combine this likelihood with the prior   
       end;
    end;
    posteriorStateKnowledge=m/sum(sum(m)); %normalize this distribution to make it a proper probability distribution.
    figure(1);mesh(posteriorStateKnowledge), axis([0 40 0 40 0 0.015]) %plot it
    figure(2);
    [a,b]=find(posteriorStateKnowledge==max(max(posteriorStateKnowledge)));  % Get the peak value; it's most likely location of the Quail.
    sest=[xLocations(a);yLocations(b)];  %A store the coordinates of this location in the bushes
    subplot(211);plot(n,sest(1),'k.');axis([0 numberOfHeardSquawks 2 4 ])
    subplot(212); plot(n,sest(2),'k.');axis([0 numberOfHeardSquawks 4 6 ])
%     pause
end;  
subplot(211); hold off;
subplot(212); hold off;