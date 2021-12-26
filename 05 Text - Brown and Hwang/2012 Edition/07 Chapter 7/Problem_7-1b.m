%This file is for Problem 7.1, Part (b)

clear all
 
Tx=[0 100;
    20 100];                    % Station coordinates
 
dt=1;                           % Sampling interval
xprior=[0; 0];                  % Initial state estimate vector
Pprior=[200^2/2 0; 0 157^2/2];   % Initial error covariance matrix
STM=[1 dt;0 1];                 % State transition matrix
 
S=400;                                  % Spectral amplitude of acceleration white noise
Q=[S*dt^3/3 S*dt^2/2; S*dt^2/2 S*dt];   % Process noise covariance Q matrix
R=4^2*eye(2);                  % Measurement noise covariance matrix
 
randn('state',13);             % Seed of randomization
store=zeros(2,21);             % Define storage
 
for t=0:20
    Rx=[200*cos(pi*t/4) 0];    % Truth position coordinates
    measnoise=randn(2,1);      % Random samples for measurement noise 
    zmeas=[norm(Rx-Tx(1,:));norm(Rx-Tx(2,:))]+4.*measnoise;  % Form noisy measurements
 
    if t==0                    % Establish initial nominal position and velocity
        Rxnom=[Rx(1)-50 0];    % based on arbitrary perturbation of initial position 
        vel=0;                 % and arbitrary zero for unknown initial velocity
    end;
    
    vtr=[Rxnom-Tx(1,:);Rxnom-Tx(2,:)];
    H=[vtr(1,1)/norm(vtr(1,:)) 0; vtr(2,1)/norm(vtr(2,:)) 0]; % Form H matrix based on nominal position
    znom=[norm(vtr(1,:)); norm(vtr(2,:))];                    % Form predicted measurements
    
    PH=Pprior*H';
    gain=PH*inv(H*PH+R);                % Compute Kalman gain
    res=zmeas-znom;                     % Subtract noisy measurement from predicted measurement to feed
                                        %    linearized model of EKF
    xupdate=xprior+gain*(res-H*xprior); % EKF state estimate update
    Rxnom(1)=Rxnom(1)+xupdate(1,1);     % Updating nominal position
    vel=vel+xupdate(2,1);               %    and velocity
    xupdate(1,1)=0;                     % Zeroing EKF state estimates after nominal updates
    xupdate(2,1)=0;
    
    IKH=eye(2)-gain*H;                  
    Pupdate=IKH*Pprior*IKH'+gain*R*gain';   % Updating error covariance matrix
    xprior=STM*xupdate;                     % Projecting EKF state estimates (this step is not really needed because of zeroed estimates)
    Pprior=STM*Pupdate*STM'+Q;              % Projecting EKF error covariance (this step is necessary)
 
    store(1,t+1)=Rx(1)-Rxnom(1);    % Save away position error
    store(2,t+1)=Rxnom(1);          % Save away estimated position
    
    Rxnom(1)=Rxnom(1)+vel;      % Projecting nominal position to next step with nominal velocity
end;
 
time=0:20;
figure(1);
plot(time,store(1,:),'-+');     % Plot position error results
axis([0 20 -60 60]);
grid;
xlabel('Time (s)');
ylabel('Error (m)');
