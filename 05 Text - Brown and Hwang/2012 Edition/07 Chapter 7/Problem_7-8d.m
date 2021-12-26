%This file is for Problem 7.8, Part (d)

clear all
 
randn('seed',12);
rand('seed',12);
N=150;              % N is the number of time steps
x=1:25000;          % The terrain profile is to be generated over 25000 m 
dt=1.0;
 
% Second-order gauss-markov model parameters for generating terrain
sigma=20.0; 
sigmasq=sigma^2;
omega0=0.01;
bsq=2*(sqrt(2))*omega0^3;
 
%Compute state model PHI and Q
F=[0 1;-omega0^2 -(sqrt(2))*omega0];
GWGT=[0 0;0 bsq*sigmasq];
A=dt*[-F GWGT;zeros(2,2) F'];
B=expm(A);
PHI=B(3:4,3:4)';
Q=PHI*B(1:2,3:4);
 
% Now generate the htrue sequences for the terrain
C=chol(Q)';
htrue=zeros(2,25000);%
htrue(1,1)=sigma*randn(1,1);
htrue(2,1)=sigma*omega0*randn(1,1);
for i=2:25000
    htrue(:,i)=PHI*htrue(:,i-1)+C*randn(2,1);
end;
h=50+htrue(1,:);
 
% Second-order gauss-markov model parameters for generating motion dynamics
sigmaz=10.0;
sigmasqz=sigmaz^2;
omega0z=0.51;
bsqz=2*(sqrt(2))*omega0z^3;
bz=sqrt(bsqz);
 
% Stack two identical models to generate horizontal and vertical components
% Compute state model PHI and Q
Fz=[0 1;-omega0z^2 -(sqrt(2))*omega0z];
Fz2=[Fz zeros(2,2); zeros(2,2) Fz];
GWGTz=[0 0;0 bsqz*sigmasqz];
GWGTz2=[GWGTz zeros(2,2); zeros(2,2) GWGTz];
Az=dt*[-Fz2 GWGTz2; zeros(4,4) Fz2'];
Bz=expm(Az);
PHIz=Bz(5:8,5:8)';
Qz=PHIz*Bz(1:4,5:8);
Cz=chol(Qz)';
 
% Generate 2nd order motion dynamics for the horizontal and vertical axes
ztrue=zeros(2,150);
ztrue(1,1)=sigmaz*randn(1,1);
ztrue(2,1)=sigmaz*omega0z*randn(1,1);
ztrue(3,1)=sigmaz*randn(1,1);
ztrue(4,1)=sigmaz*omega0z*randn(1,1);
xtrue_nom(1:2,1)=[5000; 200];
xtrue(:,1)=xtrue_nom(:,1)+[ztrue(1,1);ztrue(3,1)];
for i=2:150
    ztrue(:,i)=PHIz*ztrue(:,i-1)+Cz*randn(4,1);
    xtrue_nom(1,i)=xtrue_nom(1,i-1)+100;
    xtrue_nom(2,i)=200;
    xtrue(:,i)=xtrue_nom(:,i)+[ztrue(1,i);ztrue(3,i)];
end;
 
for k=1:N
    terrain_hgt=intpolate(x,h,xtrue(1,k));
    y(k,1)=xtrue(2,k)-terrain_hgt;  % true height above terrain
end;
yn=y+1*randn(N,1);                  % noisy height above terrain
 
% Plot terrain profile and aircraft motion profile
figure(1);
plot(x,h,'-g',xtrue(1,:),xtrue(2,:),'-b','linewidth',3);
axis([0 25000 0 300]);
xlabel('Position (m)');
ylabel('Altitude (m)');
 
%% Particle filter processing of measurements
Qmtx=5^2;                   % Defining Q process noise variance
R=1;                        % Defining R measurement noise variance
Np=100;                     % Choose 100 particles
w=ones(Np,1)/Np;            % Initialize uniform weights
xpplus=50*(97+2*rand(Np,1));    % Particle estimates uniformly-distributed between 4850 and 4950
xhat(1,:)=w(:)'*xpplus;     % Composite state estimate from weighted particle estimates
Pupdate0=sum((w.*(xpplus-xhat(1,:)))'*(xpplus-xhat(1,:))); % Composite updated covariance
for i=1:Np
    Pupdate(i)=Pupdate0;    % Assignment to each particle
end;
vv=sqrt(Qmtx).*randn(Np,1); % Generate random process noise samples for projection
xpplus=xpplus+vv;           % Project particle estimates to the next step
 
% Main recursive loop over N=150 time steps
for k=1:N    
    xprior=xpplus+101*ones(100,1);    % Projecting state estimate by nominal assumption of 101 m/s
    Pprior=Pupdate+Qmtx.*ones(1,100); % Projecting particle error covariance
    
    % Process each particle
    for i=1:Np
        Hp(i)=-(intpolate(x,h,xprior(i)+5)-intpolate(x,h,xprior(i)-5))/10;  % Linearized measurement connection
        hpred=intpolate(x,h,xprior(i));     % Calculate predicted height measurement
        yp(i)=xtrue(2,k)-hpred;             % Form measurement residual
        Pxy=Pprior(i)*Hp(i)';               % Esimate state-measurement covariance
        Pyy=(Hp(i)*Pprior(i)*Hp(i)'+R);     % Estimate measurement residual covariance
        gain(i)=Pxy/Pyy;                    % Calculate Kalman gain
        xupdate(i)=xprior(i)+gain(i)*(yn(k)-yp(i)); % Update state estimate
        Pupdate(i)=(1-gain(i)*Hp(i))^2*Pprior(i)+gain(i)^2*R;   % Update particle error covariance
        xpplus(i,1)=xupdate(i);             % Assign state estimate to particle estimate
        
        % Update particle weight
        w(i,1)=w(i,1)*exp(-0.5*(yn(k)-yp(i))*(yn(k)-yp(i))/Pyy)*exp(-0.5*(xpplus(i,1)-xprior(i))*(xpplus(i,1)-xprior(i))/Pprior(i))...
               /exp(-0.5*(xpplus(i,1)-xupdate(i))*(xpplus(i,1)-xupdate(i))/Pupdate(i)); 
    end;
    w=w/sum(w);                         % Normalize weights
    [xpplus,ww]=resample(xpplus,w,Np);  % Resample weights
    w=ww';
    xhat(k,:)=w(:)'*xpplus;             % Composite state estimate
    Pupdate0=sum((w.*(xpplus-xhat(k,:)))'*(xpplus-xhat(k,:)));  % Composite updated covariance
    for i=1:Np
        Pupdate(i)=Pupdate0;            % Assignment to each particle
    end;
    vv=sqrt(Qmtx).*randn(Np,1);         % Generate random process noise samples for projection
    xpplus=xpplus+vv;                   % Project particle estimates to the next step
end;
 
time=1:150;
figure(2);
plot(time,(xtrue(1,:)-xhat'),'x-','linewidth',2);   % Plot stored state estimates minus true values
axis([0 150 -30 30]);
grid;
xlabel('Time (s)');
ylabel('Position error (m)');
