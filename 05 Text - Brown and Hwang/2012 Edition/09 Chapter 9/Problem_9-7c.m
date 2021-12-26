%This file is for Problem 9.7, Part (c)

clear all
 
randn('seed',55);
true_range=100; % meters
c=299792458.0; % meters per second
lambda1=c/1575.42E6;  % meters
lambda2=c/1227.6E6;  % meters
elev=30*pi/180; % 30 degree SV elevation
 
xprior=zeros(3,1);
Pprior=1E6*eye(3);
Q=[1 0 0;0 0 0;0 0 0];
 
loglkhd=zeros(21,21);
store=zeros(11,2);
 
for t=1:11
    hh=cos(elev);
    meas_range=true_range*hh+1.0*randn(1); % noisy range measurement with std dev of 4m
    if t==1
        ambiguity1=round(true_range*hh/lambda1); % true value of the 1st ambiguity
        ambiguity2=round(true_range*hh/lambda2); % true value of the 2nd ambiguity
        ambiguity12 = round(true_range*hh/lambda12);
        nom_amb1=round(meas_range/lambda1);   % approx. value of the 1st ambiguity based on the noisy range measurement
        nom_amb2=round(meas_range/lambda2);   % approx. value of the 2nd ambiguity
        nom_amb12 = round(meas_range/lambda12);
    end;
    carr_noise1=0.02*randn(1);
    carr_noise2=0.02*randn(1);
    carr_phase1=true_range*hh/lambda1-ambiguity1+carr_noise1; % std dev is 5% of wavelength 1
    carr_phase2=true_range*hh/lambda2-ambiguity2+carr_noise2; % std dev is 5% of wavelength 2
 
    z=[meas_range; carr_phase1; carr_phase2]+[0; nom_amb1; nom_amb2]; % an exercise in normalizing the uncertain range of the ambiguities

end;
