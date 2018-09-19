clear
%File name is PB11_1.M and this is for Problem 11.1.

%Declare Allan variance coefficients
  h0=9.4e-20;
  h1=1.8e-19;
  h2=3.8e-21;
  time=zeros(1,61);

%Compute true Allan variance curve over the time range [0.01,1000] sec
  for t=1:61,
    tau=10^((t-21)/10);
    AllanVar(t)=h0/(2*tau)+2*h1+2/3*pi^2*h2*tau;
    time(t)=tau;
  end;

%Compute approximate Allan variance curve without flicker floor component
  for t=1:61,
    tau=10^((t-21)/10);
    AV_approx(t)=h0/(2*tau)+2/3*pi^2*h2*tau;
  end;

%Plot the true and approximate Allan variance curves for comparison
  loglog(time,AllanVar,time,AV_approx);
  axis([0.01 1000 1e-20 1e-15]);
  title ('Press ENTER to Continue');
  pause;


%Compute modified Allan variance curve h0 and h-2 parameters
  tau1=1;
  AV1=h0/(2*tau1)+2*h1+2/3*pi^2*h2*tau1;
  tau2=100;
  AV2=h0/(2*tau2)+2*h1+2/3*pi^2*h2*tau2;

  coef=[1/2    2/3*pi^2;
        1/200  200/3*pi^2];
  AV=[AV1;AV2];
  hh=inv(coef)*AV;            % solving simultaneous equations

%Compute modified Allan variance curve over [0.01 1000] sec.
  for t=1:61,
    tau=10^((t-21)/10);
    AV_modified(t)=hh(1,1)/(2*tau)+2/3*pi^2*hh(2,1)*tau;
  end;

%Plot the true and modified Allan variance curves for comparison
%Scaling below is for MATLAB 4.  Scaling was avoided for MATLAB S3.5
%because of unresolved problems when involving semilog plots.
  loglog(time,AllanVar,time,AV_modified);
  axis([0.01 1000 1e-20 1e-15]);
  title ('Press ENTER to end Problem 11.1');

%Return graphics to auto scaling mode.
  clear axis;
