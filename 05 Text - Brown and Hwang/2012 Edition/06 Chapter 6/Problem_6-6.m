%This file is for Problem 6.6

clear all
 
randn('seed',13);
dt=1;       % sampling interval
STM=1;      % state transition
Q=1;        % process noise variance
R=1.0;      % measurement noise variance
sqrtR=R^0.5;    % to scale unity noise for the measurement noise
sqrtQ=Q^0.5;    % to scale unity noise for the process noise
 
% Generating the true process, x
x=5*randn(1,1);  
for t=1:100
    x=STM*x+sqrtQ*randn(1,1);       % 1st order gauss-markov process
    H=cos(2*pi*t/200);              % time-varying measurement connection
    z(t,1)=H*x+4+sqrtR*randn(1,1);  % noisy measurement from true process
end;
 
% Kalman filter processing for a bank of 5 one-state Kalman filters
xprior=zeros(1,5);          % Initializing the bank of state estimates
Pprior=1E2*ones(1,5);       % Initializing the bank of initial error variances
loglkhd=zeros(1,5);         % Initializing the bank of log-likelihood functions
 
for t=1:100             % Iterate 100 time steps
    probtotal=0;        % Initialize probability sum variable
    H=cos(2*pi*t/200);  % time-varying measurement connection
 
    for k=1:5               % Loop to cycle through the bank of 5 filters
        PH=Pprior(1,k)*H';  
        resvar=H*PH+R;      % measurement residual variance
        gain=PH/resvar;     % Kalman gain computation      
 
        measres=(z(t,1)-k-H*xprior(1,k));               % measurement residual
        loglkhd(1,k)=loglkhd(1,k)-0.5*measres^2/resvar; % compute log-likelihood (cumulative sum of exponential part of Gaussian density). Scaling part of Gaussian density not needed because of equal variances across hypotheses
 
        xupdate(1,k)=xprior(1,k)+gain*measres;      % update state estimate
        IKH=1-gain*H;
        Pupdate=IKH*Pprior(1,k)*IKH'+gain*R*gain';  % update error variance P
 
        xprior(1,k)=STM*xupdate(1,k);       % Project ahead state estimate
        Pprior(1,k)=STM*Pupdate*STM'+Q;     % and error variance P
     
        lkhd(1,k)=exp(loglkhd(1,k));        % Form likelihood function
        probtotal=probtotal+lkhd(1,k);      % Sum likelihood function
    end;
    
    for k=1:5
        cprob(t,k)=lkhd(1,k)/probtotal;     % Compute conditional probability weights
    end;
    
end;
figure(1);
time=1:100;
plot(time,cprob(:,1),'x-',time,cprob(:,2),'o-',time,cprob(:,3),'v-',time,cprob(:,4),'+-',time,cprob(:,5),'*-');
grid;
xlabel('Steps');
ylabel('Probability');
legend('Bias 1','Bias 2','Bias 3','Bias 4','Bias 5');
