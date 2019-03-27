

%% Acknowledgements

% Written by:  StudentDave
%
% For licensing and usage questions e-mail:  scienceguy5000 at gmail. com

% Recursive Bayesian estimation example adapted from Michael A. Goodrich - September 3, 2004 - CS 470



%% Environment

close all; clear; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');

h1=figure(1); h2=figure(2);



%% Simulation Parameters

numberOfHeardSquaks=1e2;
quailLocation=[ 3; 5 ];

% The hidden quail squawks "numberOfHeardSquaks".  For each squawk the ninja hears the quail
% at some random location.  The error between where this spot is and the physical location of the
% quail will be defined as Gaussian centered at the quail's location.  The standard deviation of this
% 2-dimensional Gaussion distribution is 2 in both dimensions.

n=2*randn(2,numberOfHeardSquaks); % Create vector of iterative squawks with a standard deviation of 2.
x=zeros(2,numberOfHeardSquaks); % x will be the initialized variable for where the ninja thinks he hears the quail.

%center the gaussian random swarks around the point where the quail
%actually is. and plot
figure(1);

%make the plot prettier
h=plot(quailLocation(1),quailLocation(2),'r.');  % Plot where the quail actually is
set(h,'markersize',40, 'linewidth',3); % make pretty
axis([0,10,0,10]);  % make pretty
hold off;
hold on
for i=1:numberOfHeardSquaks
    x(:,i)=quailLocation+n(:,i);
    plot(x(1,i),x(2,i),'k.','markersize',10);
end;




%now, we do Ninja bayes!

%define the locations the quail can be at
Sa=[2:0.05:4];
Sb=[4:0.05:6];


%no bias, no prior knowledge, uniform distribution
L=length(Sa);
Pr=ones(L,L); % Initialize all one --> uniform prior
Po=ones(L,L); %duplicate for iterative update

%bias the prior toward the quail
%figure(1);mesh(centered_prior), axis([0 40 0 40 0 0.015])
%Po = centered_prior;
%Pr = centered_prior;

%bias the prior away from the quail
%figure(1);mesh(off_centered_prior), axis([0 40 0 40 0 0.015])
%Po = off_centered_prior;
%Pr = off_centered_prior;

Pr=Pr/sum(sum(Pr)); % Turn the prior into a pmf by dividing by the sum.
Po=Po/sum(sum(Po)); % Each value is now 1/(number of states), so it all sums to one.
figure(1);clf;mesh(Po), axis([0 40 0 40 0 0.015])






%%iterative bayes

[a,b]=find(Po==max(max(Po)));  % Pull out the indices at which Po achieves its max to start.
sest=[Sa(a);Sb(b)];  % The best estimate of the true state to start.
figure(1);
clf
figure(2);
clf
subplot(211); plot(1,sest(1)); hold on;
line([1,numberOfHeardSquaks],[quailLocation(1),quailLocation(1)]); % Draw a line at the location of the x dimension.
subplot(212); plot(1,sest(2)); hold on;
line([1,numberOfHeardSquaks],[quailLocation(2),quailLocation(2)]); % Draw a line at the location of the y dimension.

K=[4,0;0,4]; % covariance matrix for making a 2-D gaussian
for (n=2:length(x));
    Pr=Po; %store the posterior to the prior.
    m=0*Pr;   
    %likelihood
    % look at each location, assume that the given location is
    % is there quail is, and get the likelihood of the data x(:,n) assuming
    % 2-d gaussian noise
    for (i=1:length(Pr))
       for (j=1:length(Pr))
           me=[Sa(i);Sb(j)];
           m(i,j) = 1/sqrt((2*pi)^2*det(K)) * exp(-(x(:,n)-me)'*inv(K)*(x(:,n)-me)/2); %Compute likelihood           
           m(i,j) = m(i,j) * Pr(i,j); % Combine this likelihood with the prior   
       end;
    end;
    Po=m/sum(sum(m)); %normalize this distribution to make it a proper probability distribution.
    figure(1);mesh(Po), axis([0 40 0 40 0 0.015]) %plot it
    figure(2);
    [a,b]=find(Po==max(max(Po)));  % Get the peak value; it's most likely location of the Quail.
    sest=[Sa(a);Sb(b)];  %A store the coordinates of this location in the bushes
    subplot(211);plot(n,sest(1),'k.');axis([0 numberOfHeardSquaks 2 4 ])
    subplot(212); plot(n,sest(2),'k.');axis([0 numberOfHeardSquaks 4 6 ])
%     pause
end;  
subplot(211); hold off;
subplot(212); hold off;