% This file is for Problem 9.1, Part (b)
% 1-D positioning problem - Problem 9.1
% Monte Carlo simulation of position estimation of a wandering location
 
clear all
 
randn('seed',51);             %fix the seed of randomization
Nom_pos=[0,6378000];
SV_radius=26560000;
SV_init_elev=[30*pi/180; 135*pi/180];
SV_elev_rate=[2*pi/43200; -2*pi/43200];
 
True_pos(1,:)=Nom_pos+[100*randn(1),0]; % Set initial position
clock_err=10000*randn(1);               % Set initial clock error
for t=0:10000
    SV_elev=SV_init_elev+SV_elev_rate*t;
    SV_pos=[SV_radius.*cos(SV_elev) SV_radius.*sin(SV_elev)]; %rows are different SVs; columns are x-y components
    drnvtr=ones(2,1)*True_pos(t+1,:)-SV_pos;
    pseudorange(t+1,:)=[norm(drnvtr(1,:))+clock_err; norm(drnvtr(2,:))+clock_err]+10*randn(2,1);
    True_pos(t+2,:)=True_pos(t+1,:)+[2*randn(1) 0]; % Project position for next cycle (random walk model)
    clock_err=clock_err+0.1*randn(1);               % Project clock error for next cycle (random walk model)
end;
 
Q=[2^2 0;0 0.1^2];
R=10^2*eye(2);
xprior=[0; 0];
Pprior=[100^2 0;0 10000^2];
for t=0:10000
    SV_elev=SV_init_elev+SV_elev_rate*t;
    SV_pos=[SV_radius.*cos(SV_elev) SV_radius.*sin(SV_elev)]; %rows are different SVs; columns are x-y components
    
    drnvtr=ones(2,1)*Nom_pos-SV_pos;
    H=[drnvtr(1,1)/norm(drnvtr(1,:)) 1; drnvtr(2,1)/norm(drnvtr(2,:)) 1];
    PH=Pprior*H';
    gain=PH/(H*PH+R);                                   %compute Kalman gain
    pred_range=[norm(drnvtr(1,:)); norm(drnvtr(2,:))];  %compute predicted pseudorange
    meas_res=pseudorange(t+1,:)'-pred_range-H*xprior;   %form measurement residual
    xupdate=xprior+gain*meas_res;                       %update state estimate
    Pupdate=(eye(2)-gain*H)*Pprior;                     %update P matrix
    xprior=xupdate;                                     %project state estimate
    Pprior=Pupdate+Q;                                   %project P matrix
    store(t+1,1)=xupdate(1,1)+Nom_pos(1,1)-True_pos(t+1,1); %save position error
    store(t+1,2)=Pupdate(1,1)^0.5;                          %save (1,1) term of Pmatrix
    store(t+1,3)=xupdate(1,1);                              %save position estimate (deviation from nominal)
    store(t+1,4)=True_pos(t+1,1);                           %save true position (deviation from nominal)
end;
figure(1);
time=0:10000;
plot(time,store(:,1),time,store(:,2),'r--',time,-store(:,2),'r--','linewidth',2);
axis([0 10000 -30 30]);
xlabel('Time (s)');
ylabel('Position error (m)');
