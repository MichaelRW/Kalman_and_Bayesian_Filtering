%This file is for Problem 6.5

clear all
 
randn('seed',3);
dt=1;               % Sampling interval
STM=[1 dt; 0 1];    % State transition matrix
S1=0.01;            % Spectral amplitude of process noise for Hypothesis 1
S2=1;               % Spectral amplitude of process noise for Hypothesis 2
Q1=[S1*dt^3/2 S1*dt^2/2; S1*dt^2/2 S1*dt];  % Process noise covariance for Hypo. 1
Q2=[S2*dt^3/2 S2*dt^2/2; S2*dt^2/2 S2*dt];  % Process noise covariance for Hypo. 2
H=[1 0];            % Measurement connection
R=1.0;              % Measurement noise variance
sqrtR=R^0.5;        % for scaling unity noise to get measurement noise
sqrtQ2=chol(Q2)';   % for premultiplying 2-tuple unit noise vector to get process noise
 
% Generating the true process, x (2-tuple)
x(:,1)=zeros(2,1);                      % Setting up initial states
z(1,1)=H*x(:,1)+sqrtR*randn(1,1);       % Noisy measurement related to initial states
for t=2:20
    x(:,t)=STM*x(:,t-1)+sqrtQ2*randn(2,1);  % Integrated random walk process (Q2 used means Hypothesis 2 should be the correct one)
    z(t,1)=H*x(:,t)+sqrtR*randn(1,1);       % Noisy measurement of true process
end;
 
% Kalman filter processing- setting up two 2-tuple Kalman filter hypotheses
xprior1=zeros(2,1);     % Initializing state and error covariance for Hypo. 1
Pprior1=1E6*eye(2);
xprior2=zeros(2,1);     % Initializing state and error covariance for Hypo. 2
Pprior2=1E6*eye(2);
loglkhd1=0;             % Initializing log-likelihood functions for hypotheses
loglkhd2=0;
sf1=1;
sf2=1;
 
for t=1:20              % Iterate over 20 time steps
    PH1=Pprior1*H';     % Processing Kalman filter hypothesis 1
    resvar1=H*PH1+R;                % measurement residual variance
    gain1=PH1/resvar1;              % Compute Kalman gain
    measres1=(z(t,1)-H*xprior1);    % measurement residual
    loglkhd1=loglkhd1-0.5*measres1^2/resvar1;   % log-likelihood for Hypo. 1 (cumulative sum of exponential part of Gaussian density)
    sf1=sf1/((2*pi*resvar1)^0.5);	    % cumulative product of scaling part of Gaussian density for Hypo. 1
    xupdate1=xprior1+gain1*measres1;    % update state estimate vector
    IKH1=eye(2)-gain1*H;
    Pupdate1=IKH1*Pprior1;          % update error covariance P
    xprior1=STM*xupdate1;           % Project ahead the state estimate vector
    Pprior1=STM*Pupdate1*STM'+Q1;   % Project ahead the error covariance P
    
    PH2=Pprior2*H';     % Same processing steps as above for Kalman filter Hypo. 2
    resvar2=H*PH2+R;
    gain2=PH2/resvar2;
    measres2=(z(t,1)-H*xprior2);
    loglkhd2=loglkhd2-0.5*measres2^2/resvar2;
    sf2=sf2/((2*pi*resvar2)^0.5);
    xupdate2=xprior2+gain2*measres2;
    IKH2=eye(2)-gain2*H;
    Pupdate2=IKH2*Pprior2;
    xprior2=STM*xupdate2;
    Pprior2=STM*Pupdate2*STM'+Q2;
    
    lkhd1=sf1*exp(loglkhd1);    % Form likelihood function for Hypo. 1
    lkhd2=sf2*exp(loglkhd2);    % Form likelihood function for Hypo. 2 
    cprob1(t,1)=lkhd1/(lkhd1+lkhd2);    % Compute conditional probability weight
    cprob2(t,1)=lkhd2/(lkhd1+lkhd2);    % Compute conditional probability weight
    
    xout(t,1)=cprob1(t,1)*xupdate1(1,1)+cprob2(t,1)*xupdate2(1,1);  % Form weighted state estimate
end;
figure(1);
time=1:20;
plot(time,cprob1,'x-',time,cprob2,'o-');    % Plot conditional probability profile for both hypotheses
grid;
xlabel('Steps');
ylabel('Probability');
legend('Hypothesis 1: S=0.01','Hypothesis 2: S=1.0','Location','SouthEast');
 
figure(2);
time=1:20;
plot(time,xout(:,1),'x-',time,x(1,:),'o-'); % Plot weighted state estimate alongside true process (first state only)
grid;
axis([0 20 -10 70]);
xlabel('Steps');
ylabel('Probability');
legend('Weighted estimate of State 1','True value of State 1','Location','NorthWest');
