clear
% File name is PB11_5.M and this is for Problem 11.5.

% First, specify process and measurement noise parameters.
  phi=1;
  q=100;
  sigma_process=sqrt(q)
  r=225;
  sigma_meas=sqrt(r);
  pminus=100;
  TotalSteps=300;
  time=0:1:TotalSteps;

% Generate random process and noisy measurements with a fixed seed 0 (may
% be changed if desired).
  rand ('seed',0);
  rand ('normal');
  xtrue=zeros(1,TotalSteps+1);
  z=zeros(1,TotalSteps+1);
  xtrue(1)=sqrt(pminus)*rand;

  disp ('Running...');
  for k=1:TotalSteps+1
    time(k)=k-1;
    if k==1,
      xtrue(k)=sigma_process*rand;
    else
      xtrue(k)=xtrue(k-1)+sigma_process*rand;
    end;
    angle=75+(k-1)*3600/43082;
    h=-cos(angle*pi/180);
    z(k)=h*xtrue(k)+sigma_meas*rand;
  end;

% Process measurements generated above using Kalman filter and
% compute its estimation error.
  xminus=0;
  xplus=zeros(1,TotalSteps+1);
  for k=1:TotalSteps+1
    angle=75+(k-1)*3600/43082;
    h=-cos(angle*pi/180);
    gain=pminus*h/(h*pminus*h+r);
    xplus(k)=xminus+gain*(z(k)-h*xminus);
    pplus=(1-gain*h)*pminus;
    xminus=phi*xplus(k);
    pminus=phi*pplus*phi+q;
  end;
  KF_error=xplus-xtrue;     %generate Kalman filter estimation error

% Process measurements generated above using 10-sec Moving Average
% filter and compute the estimation error.
  mov_avg=zeros(1,TotalSteps+1);
  sum=0;
  for k=1:TotalSteps+1
    sum=sum+z(k);
    if k<=10,
      z(k)=sum/k;
    else
      z(k)=sum/10;
      sum=sum-z(k-9);
    end;
    angle=75+(k-1)*3600/43082;
    h=-cos(angle*pi/180);
    mov_avg(k)=z(k)/h;
  end;
  MA_error=mov_avg-xtrue;    %generate moving average estimation error

% Plot the estimation errors obtained from the two different ways
% of filtering the same data.
  axis ([0 TotalSteps -500 500]);
  plot(time,KF_error,time,MA_error);
  axis ([0 TotalSteps -500 500]);

%Return graphics to auto scaling mode.
  clear axis;
