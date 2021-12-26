%This file is for Problem 9.8, Part (b)

clear all;
randn('seed',1);       % Set a fixed seed of randomization
dt = 1;                % Sampling interval of 1 second
STMone = [1 dt; 0 1];  % State transition submatrix
STM = [STMone zeros(2,4); zeros(2,2) STMone zeros(2,2); zeros(2,4) STMone]; % Assembling full state transition matrix
S = 2;                                       % Spectral amplitude for motion dynamics
Qrange = [S*dt^3/3 S*dt^2/2; S*dt^2/2 S*dt]; % Q matrix for motion dynamics
Sf_TCXO = 1.348E-4;          % Spectral amplitude for white noises in TCXO model
Sg_TCXO = 3.548E-6;
QclkTCXO = [Sf_TCXO*dt+Sg_TCXO*dt^3/3  Sg_TCXO*dt^2/2; Sg_TCXO*dt^2/2  Sg_TCXO*dt]; % Q matrix for TCXO model
Sf_Rb = 9E-6;                % Spectral amplitude for white noises in Rb model
Sg_Rb = 1.8E-12;
QclkRb = [Sf_Rb*dt+Sg_Rb*dt^3/3  Sg_Rb*dt^2/2; Sg_Rb*dt^2/2  Sg_Rb*dt]; % Q matrix for Rb model
Q = [Qrange zeros(2,4); zeros(2,2) QclkTCXO zeros(2,2); zeros(2,4) QclkRb]; % Assembling full Q matrix
H1=[1 0 1 0 -1 0];    % H matrix for aircraft-to-ground pseudorange measurement
H2=[1 0 -1 0 1 0];    % H matrix for ground-to-aircraft pseudorange measurement
R=2^2;                % R matrix for pseudorange measurement noise covariance
 
sqrtQ_TCXO=chol(QclkTCXO)';  % Scale factor for generating TCXO process noise
sqrtQ_Rb=chol(QclkRb)';      % Scale factor for generating Rb process noise
gnd_x=-10000;    % Ground station x-coordinate
gnd_y=0;         % Ground station y-coordinate
rw_x=0;      % Initial value of x-position random walk
rw_y=0;      % Initial value of y-position random walk
xprior=zeros(6,1); % Initial KF state estimate
xprior(2,1)=70;    % Setting an approximate initial range velocity state (State 2)
Pprior=1E6*eye(6); % Setting the initial error covariance P
clk_err_TCXO=[1E3*randn(1,1); 10*randn(1,1)]; % Initial TCXO states
clk_err_Rb=[1E3*randn(1,1); 0];               % Initial Rb states
for t=1:200
    nom_x=0;                  % Nominal position changes in x and y
    nom_y=-10000+100*t;
    rw_x=rw_x+1*randn(1,1);   % Random walk perturbation to nominal position in x and y
    rw_y=rw_y+1*randn(1,1);
    pos_x=nom_x+rw_x;         % True position of aircraft in x and y
    pos_y=nom_y+rw_y;
    range=norm([pos_x pos_y]-[gnd_x gnd_y]); % True range between ground and aircraft
    nom_range=norm([nom_x nom_y]-[gnd_x gnd_y]); % Nominal range between ground and nominal aircraft position
    clk_err_TCXO=STMone*clk_err_TCXO+sqrtQ_TCXO*randn(2,1);  % Generating TCXO random process sequence
    clk_err_Rb=STMone*clk_err_Rb+sqrtQ_Rb*randn(2,1);        % Generating Rb random process sequence
    
    if mod(t,10) == 1 % Time for aircraft-to-ground pseudorange measurement update
        meas=range+clk_err_TCXO(1,1)-clk_err_Rb(1,1)+2*randn(1,1); % Aircraft-to-ground pseudorange measurement
        PH=Pprior*H1';
        gain=PH/(H1*PH+R);                          % KF gain
        xupdate=xprior+gain*(meas-H1*xprior);       % KF state estimate update
        I_KH=eye(6)-gain*H1;
        Pupdate=I_KH*Pprior*I_KH'+gain*R*gain';     % KF error covariance update
    else
        if mod(t,10) == 5 % Time for ground-to-aircraft pseudorange measurement update
            meas=range-clk_err_TCXO(1,1)+clk_err_Rb(1,1)+2*randn(1,1); % Ground-to-aircraft pseudorange measurement
            PH=Pprior*H2';
            gain=PH/(H2*PH+R);                      % KF gain
            xupdate=xprior+gain*(meas-H2*xprior);   % KF state estimate update
            I_KH=eye(6)-gain*H2;
            Pupdate=I_KH*Pprior*I_KH'+gain*R*gain'; % KF error covariance update
        else
            xupdate=xprior;                         % Trivial no-measurement state update otherwise
            Pupdate=Pprior;                         % Trivial no-measurement covariance update otherwise
        end;
    end;
    store(t,1)=xupdate(1,1)-range;                                            % Saving range error 
    store(t,2)=xupdate(3,1)-xupdate(5,1)-(clk_err_TCXO(1,1)-clk_err_Rb(1,1)); % Saving relative timing error
    store(t,3)=Pupdate(1,1)^0.5;                                              % Saving range standard deviation 
    store(t,4)=(Pupdate(3,3)+Pupdate(5,5)-2*Pupdate(3,5))^0.5;                % Saving relative timing error standard dev.
    xprior=STM*xupdate;                     % State estimate projection to the next cycle
    Pprior=STM*Pupdate*STM'+Q;              % Error covariance projection to the next cycle
end;
 
index=0:199;
plot(index,store(:,2),'-',index,store(:,4),'r--',index,-store(:,4),'r--','LineWidth',2); % Generate plots of timing error and timing error standard dev.
xlabel('Time (s)');
ylabel('Relative Timing Error (m)');
axis([1 200 -10 10]);
grid;
