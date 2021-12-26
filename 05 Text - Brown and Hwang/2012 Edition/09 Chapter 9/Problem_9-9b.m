%This file is for Problem 9.9, Part (b)

clear all
 
dt=1;   % Sampling interval is one second
 
% Ha matrix consists of unit direction vectors from Location A to the 8 satellites
Ha=[-0.2031    0.6498   -0.7325   -1.0000 0;
    -0.5305    0.0999   -0.8418   -1.0000 0;
     0.6336   -0.7683   -0.0913   -1.0000 0;
     0.4705    0.0745   -0.8792   -1.0000 0;
     0.0955    0.6480   -0.7557   -1.0000 0;
     0.7086   -0.4356   -0.5551   -1.0000 0;
    -0.6270   -0.4484   -0.6371   -1.0000 0;
    -0.8651   -0.4877   -0.1170   -1.0000 0];  
 
% Hb matrix is a submatrix of Ha from Location B to 4 of the 8 satellites
Hb=[Ha(1,:);Ha(2,:);Ha(4,:);Ha(5,:)];
 
 
Pprior=1E4*eye(10); % Specify initial P matrix
Rpr=10;             % Set the measurement error R matrix
 
Sf=1.348E-4;  % Spectral amplitudes of white noises for receiver clock
Sg=3.548E-6;
Qclk=[Sf*dt+Sg*dt^3/3 Sg*dt^2/2; Sg*dt^2/2 Sg*dt];  % Q submatrix for clock dynamics
Qobs=10*eye(5);     % Q submatrix for position dynamics
Qobs(4:5,4:5)=Qclk; % Qobs is the full Q submatrix for each location
Q=[Qobs zeros(5,5); zeros(5,5) Qobs];   % Full Q matrix for two locations
 
STMobs=eye(5);
STMobs(4,5)=dt;   % State transition submatrix for each location
STM=[STMobs zeros(5,5); zeros(5,5) STMobs]; % Full state transition matrix for two locations
 
R=Rpr*eye(12); % R matrix of GPS pseudorange measurements: 8 SVs at A and 4 SVs at B
 
Rclk=4;  % Measurement variance Rclk of the relative timing measurement
 
H=[Ha zeros(8,5); zeros(4,5) Hb];  % Full H matrix with 8 rows for A and 4 rows for B
Hclk=[0 0 0 1 0 0 0 0 -1 0];  % Additional Hclk with 1 row for the relative timing measurement
    
for t=1:120
    PH=Pprior*H';      % KF gain computations
    gain=PH/(H*PH+R);
    Pupdate=(eye(10)-gain*H)*Pprior;  % KF covariance update
    
    if t > 100              % Processing relative timing measurement after t = 100
        Pprior=Pupdate;     % Supplemental processing -- trivial same-cycle projection of updated P
        PH=Pprior*Hclk';    % KF gain computation for supplemental measurement
        gain=PH/(Hclk*PH+Rclk);
        Pupdate=(eye(10)-gain*Hclk)*Pprior;  % KF covariance update
    end;
    
    store(t,1)=sqrt(Pupdate(1,1)+Pupdate(2,2)+Pupdate(3,3));  % Compute and save 3D position error for A
    store(t,2)=sqrt(Pupdate(6,6)+Pupdate(7,7)+Pupdate(8,8));  % Compute and save 3D position error for B
    
    Pprior=STM*Pupdate*STM'+Q;  % Project P to the next cycle
end;
 
time=1:120;
plot(time,store(:,1),'-x',time,store(:,2),'-x');  % Plot results
axis([0 120 0 100]);
grid;
xlabel('Time(sec)');
ylabel('3D position error(m)');
legend('Location A','Location B');
