%This file is for Problem 9.1, Part (a)
% 1-D positioning problem 
% Error covariance analysis of a nominal location
 
clear all
 
Nom_pos=[0,6378000];
SV_radius=26560000;
SV_init_elev=[30*pi/180; 135*pi/180];
SV_elev_rate=[2*pi/43200; -2*pi/43200];
 
Q=[2^2 0;0 0.1^2];
R=10*eye(2);
Pprior=10000*eye(2);
for t=0:10000
    SV_elev=SV_init_elev+SV_elev_rate*t;
    SV_pos=[SV_radius.*cos(SV_elev) SV_radius.*sin(SV_elev)]; 
%rows are different SVs; columns are x-y components
    
    drnvtr=ones(2,1)*Nom_pos-SV_pos;
    H=[drnvtr(1,1)/norm(drnvtr(1,:)) 1; drnvtr(2,1)/norm(drnvtr(2,:)) 1];
    PH=Pprior*H';
    gain=PH/(H*PH+R);               %compute Kalman gain
    Pupdate=(eye(2)-gain*H)*Pprior; %update P matrix
    Pprior=Pupdate+Q;               %project P matrix
    
    store(t+1,1)=Pupdate(1,1)^0.5;  %save away (1,1) term of P matrix
end;
figure(1);
time=0:10000;
plot(time,store,'linewidth',2);
xlabel('Time (s)');
ylabel('Position error (m)');
