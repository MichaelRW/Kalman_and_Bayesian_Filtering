clear
% File name is PB11_6.M and this is for Problem 11.6.

% First, specify the filter parameters.

  c=3E8;
  dt=1;
  TotalSteps=60;
  state_dim=5;
  PHI=eye(state_dim);
  PHI(4,5)=dt;

  H=[-0.3523 -0.0495  0.9346  1  0;
      0.6199  0.7406  0.2593  1  0;
     -0.9506 -0.2553  0.1764  1  0;
      0.9613  0.2129  0.1747  1  0];

  R=225*eye(4);

% Set up Q matrix with values given for part (a).
  Sp=1.0;
  Sf=0.4E-18*c^2;
  Sg=1.58E-18*c^2;
  Q=[Sp*dt  0       0         0               0;
     0    Sp*dt     0         0               0;
     0      0     Sp*dt       0               0;
     0      0       0  Sf*dt+(Sg*dt^3)/3  (Sg*dt^2)/2;
     0      0       0      (Sg*dt^2)/2      Sg*dt];

% Set up a loop to solve for the three parts of the problem.
  for part=1:3
    if part==2,
      Sp=10;
      Q(1,1)=Sp*dt;
      Q(2,2)=Sp*dt;
      Q(3,3)=Sp*dt;
    else
      if part==3,
        R=2250*eye(4);
      end
    end;

    PMINUS=[10000 0 0 0 0;
            0 10000 0 0 0;
            0 0 10000 0 0;
            0 0 0 90000 0;
            0 0 0  0  900];

    time=0:1:TotalSteps;
    Save_P=zeros(3,TotalSteps+1);

    disp ('Running...');

% Main recursive loop
    for k=1:TotalSteps+1
      time(k)=k-1;
      GAIN=PMINUS*H'*inv(H*PMINUS*H'+R);
      PPLUS=(eye(state_dim)-GAIN*H)*PMINUS;
      PPLUS=(PPLUS+PPLUS')/2;
      PMINUS=PHI*PPLUS*PHI'+Q;
      PMINUS=(PMINUS+PMINUS')/2;
      Save_P(1,k)=PPLUS(1,1);
      Save_P(2,k)=PPLUS(2,2);
      Save_P(3,k)=PPLUS(3,3);
    end;

    % Plot the east, north, vertical position components in Save_P.
    axis([0 TotalSteps 0 1000]);
    plot(time,Save_P(1,:),time,Save_P(2,:),time,Save_P(3,:));
    axis([0 TotalSteps 0 1000]);
    if part==3,
      title ('Press ENTER to end Problem 11.6');
    else
      title ('Press ENTER to Continue');
    end;
    pause;
  end;

%Return graphics to auto scaling mode.
  clear axis;
