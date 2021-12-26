%This file is for Problem 9.7, Part (d)

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
    H=[hh 0 0; hh/lambda1 1 0; hh/lambda2 0 1];
    R=[1^2 0 0; 0 0.02^2 0; 0 0 0.02^2];
    
    gain=Pprior*H'*inv(H*Pprior*H'+R);         % Kalman filter estimating the range value as well as the ambiguity pair
    xupdate=xprior+gain*(z-H*xprior);          % Updating the state estimate vector
    I_KH=eye(3)-gain*H;
    Pupdate=I_KH*Pprior*I_KH'+gain*R*gain';    % Updating the error covariance matrix
    
    xprior=xupdate;                            % Projecting the state estimate vector
    Pprior=Pupdate+Q;                          % Projecting the updated error covariance (states are constant biases)
 
    crossterm=inv(Pupdate(2:3,2:3))*xupdate(2:3,1);      % Reconstructing the likelihood weighting matrix 
    W=[inv(Pupdate(2:3,2:3)) crossterm; crossterm' 0];   % for the Magill adaptive filter
    max=-1e10;
    for u1=-10:1:10
        for u2=-10:1:10
            loglkhd(u1+11,u2+11)=-0.5*[u1 u2 1]*W*[u1 u2 1]';
            if loglkhd(u1+11,u2+11) > max
                max=loglkhd(u1+11,u2+11);
            end;
        end;
    end;
    sum1=0;    
    for u1=-10:1:10
        for u2=-10:1:10
            sum1=sum1+exp(loglkhd(u1+11,u2+11)-max);
        end;
    end;
    wt=zeros(21,21);
    for u1=-10:1:10
        for u2=-10:1:10
            wt(u1+11,u2+11)=exp(loglkhd(u1+11,u2+11)-max)/sum1;
        end;
    end;
    
    prob=exp(loglkhd(ambiguity1-nom_amb1+11,ambiguity2-nom_amb2+11)-max)/sum1;  % Isolating the probability associated with the "correct hypothesis"
    store(t,1)=prob;
 
    xpos=zeros(21,21);
    Hc=[0 1 0; 0 0 1];                                                     % Isolating the state estimate of the true range value after constraining
    gain_c=Pprior*Hc'*inv(Hc*Pprior*Hc');                                  % the integer estimate to be the hypothesized ambiguity pair
    for u1=-10:1:10
        for u2=-10:1:10
            xupdate_c=xprior+gain_c*(-[u1;u2]-Hc*xprior); 
            xpos(u1+11,u2+11)=xupdate_c(1,1);                                     % The first term of the updated constrained state estimate vector is the
                                                                           % estimate of the true range value associated with this hypothesis
        end;
    end;
    store(t,2)=sum(sum(wt.*xpos));
                                                                           
end;
figure(1);
index=0:1:10;
plot(index,store(:,1),'-+','LineWidth',2);                                                % Plot the time profile of the probability associated with the
axis([0 10 0 1]);                                                          % correct hypothesis (ambiguity pair) for the 2D hypothesis space model
grid;
xlabel('Time steps');
ylabel('Probability');
figure(2);
plot(index,store(:,2),'-x','LineWidth',2);                                                % Plot the time profile of the first term of the updated constrained
axis([0 10 99.9 100.1]);                                                   % state estimate vector associated with the correct hypothesis (amb. pair)
grid;
xlabel('Time steps');
ylabel('Baseline solution (m)');
