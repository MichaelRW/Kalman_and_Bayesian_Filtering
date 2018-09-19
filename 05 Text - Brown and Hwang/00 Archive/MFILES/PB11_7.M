clear
% File name is PB11_7.M and this is for Problem 11.7.

% First, specify the filter parameters.

  c=3E8;    % speed of light
  dt=1;     % time increment
  TotalSteps=200;

% Part (a): Model using differenced measurements.
  state_dim=3;
  PHI3=eye(state_dim);
  Sp=1.0;
  Q3=[Sp*dt  0       0;
      0    Sp*dt     0;
      0      0     Sp*dt];

  H5=[-0.3523 -0.0495  0.9346  1  0;
       0.6199  0.7406  0.2593  1  0;
      -0.9506 -0.2553  0.1764  1  0;
       0.9613  0.2129  0.1747  1  0];
  H3=zeros(3);
  H3(1,:)=H5(1,1:3)-H5(2,1:3);
  H3(2,:)=H5(2,1:3)-H5(3,1:3);
  H3(3,:)=H5(3,1:3)-H5(4,1:3);

  R3=225*2*eye(3);
  R3(1,2)=-225;
  R3(2,3)=-225;
  R3(2,1)=-225;
  R3(3,2)=-225;

  PMINUS3=[10000 0 0;
           0 10000 0;
           0 0 10000];

  time=0:1:TotalSteps;
  Save_P=zeros(2,TotalSteps+1);

  disp ('Running Part(a)...');

% Main recursive loop.
  for k=1:TotalSteps+1
    time(k)=k-1;
    GAIN3=PMINUS3*H3'*inv(H3*PMINUS3*H3'+R3);
    PPLUS3=(eye(state_dim)-GAIN3*H3)*PMINUS3;
    PPLUS3=(PPLUS3+PPLUS3')/2;
    PMINUS3=PHI3*PPLUS3*PHI3'+Q3;
    PMINUS3=(PMINUS3+PMINUS3')/2;
    Save_P(1,k)=PPLUS3(1,1)+PPLUS3(2,2)+PPLUS3(3,3);
  end;

% Part (b): Model using non-differenced measurements;
%           same as in Problem 11.6(a).
  state_dim=5;
  PHI5=eye(state_dim);
  Sp=1.0;
  Sf=0.4E-18*c^2;
  Sg=1.58E-18*c^2;
  Q5=[Sp*dt  0       0         0               0;
      0    Sp*dt     0         0               0;
      0      0     Sp*dt       0               0;
      0      0       0  Sf*dt+(Sg*dt^3)/3  (Sg*dt^2)/2;
      0      0       0      (Sg*dt^2)/2      Sg*dt];

  PMINUS5=[10000 0 0 0 0;
           0 10000 0 0 0;
           0 0 10000 0 0;
           0 0 0 90000 0;
           0 0 0  0  900];

  R5=225*eye(4);

  disp ('Running Part(b)...');

% Main recursive loop.

  for k=1:TotalSteps+1
    time(k)=k-1;
    GAIN5=PMINUS5*H5'*inv(H5*PMINUS5*H5'+R5);
    PPLUS5=(eye(state_dim)-GAIN5*H5)*PMINUS5;
    PPLUS5=(PPLUS5+PPLUS5')/2;
    PMINUS5=PHI5*PPLUS5*PHI5'+Q5;
    PMINUS5=(PMINUS5+PMINUS5')/2;
    Save_P(2,k)=PPLUS5(1,1)+PPLUS5(2,2)+PPLUS5(3,3);
  end;
fig1 = figure;
% Compare plots between models used in parts (a) and (b).
  axis([0 200 0 100]);
  plot(time,Save_P(1,:),time,Save_P(2,:));
  axis([0 200 0 100]);

%Return graphics to auto scaling mode.
  clear axis;
pause;
close (fig1);
