clear
%File name is PB11_2B.M and this is for Problem 11.2 part (b).

%The purpose of this m-file is to compute the optimal Kalman filter
%error covariance over m steps.  (Note that in the notation used
%here m does not count the measurement update at k=0, i.e., there
%are a total of m+1 measurements to process.)

%The INS phi and q matrices will be obtained from m-file insphiq.m.
%This is a special file for the 9-state model given by Eq.(10.2.10).
%The variable names produced by insphiq.m are PHI and Q.

%Specify global variables and key parameters needed for phi/q calculation.

	global PHI Q dt omegax omegay wacc wgyro

	dt=1;
	omegax=-.0000157;
	omegay=.0000727;
	wacc=.0036;
	wgyro=2.35e-11;
	c=3e8;


%The measurement sequence starts at k=0 and ends at k=m.

	m=1000;
	disp ('How many steps (default is 1000)');
	minput = input ('? ');
	if minput >= 0,
  	  m = minput;
	end;
	

%Now set up row matrices with the correct dimensions for your outputs.

	time=0:1:m;


%Start by specifying one of the three different sub-problems
%illustrated by this example.

        disp ('Running...');


%PHI11 and Q11 matrices are made up of 9x9-submatrices phi and q
%obtained from the script file, insphiq.

	insphiq;
	nullvtr=zeros(9,1);
	PHI11=[PHI  nullvtr  nullvtr;
	       nullvtr'   1     dt;
	       nullvtr'   0     1];
	
	Sf=0.4E-18*c^2;
	Sg=1.58E-18*c^2;
	Q11=[Q  nullvtr  nullvtr;
	     nullvtr'  Sf*dt+(Sg*dt^3)/3  (Sg*dt^2)/2;
	     nullvtr'     (Sg*dt^2)/2      Sg*dt];


%Specify the initial P matrix as PMINUS.

        PMINUS= [1E6 0  0  0  0  0  0  0  0  0  0;
                 0  1E2 0  0  0  0  0  0  0  0  0;
                 0  0  6.4E-3 0  0  0  0  0  0  0  0;
                 0  0  0  1E6 0  0  0  0  0  0  0;
                 0  0  0  0  1E2 0  0  0  0  0  0;
                 0  0  0  0  0  6.4E-3 0  0  0  0  0;
                 0  0  0  0  0  0  1E6 0  0  0  0;
                 0  0  0  0  0  0  0  1E2 0  0  0;
                 0  0  0  0  0  0  0  0  6.4E-3 0  0;
		 0  0  0  0  0  0  0  0  0  1E6 0;
		 0  0  0  0  0  0  0  0  0  0  1E4];
     	state_dim = 11;


%Specify the numerical values for H and R.

	H=[-0.3523 0  0 -0.0495 0  0  0.9346 0  0  1  0;
	    0.6199 0  0  0.7406 0  0  0.2593 0  0  1  0;
	   -0.9506 0  0 -0.2553 0  0  0.1764 0  0  1  0;
	    0.9613 0  0  0.2129 0  0  0.1747 0  0  1  0;
	      0    0  1    0    0  0    0    0  0  0  0;
	      0    0  0    0    0  1    0    0  0  0  0;
	      0    0  0    0    0  0    0    0  1  0  0];


	R=[225 0  0  0  0   0   0; 
	    0 225 0  0  0   0   0;
	    0  0 225 0  0   0   0;
	    0  0  0 225 0   0   0;
	    0  0  0  0 1e-4 0   0;
	    0  0  0  0  0 1e-4  0;
	    0  0  0  0  0   0 1e-4];


%Now set up row matrices with the correct dimensions for the outputs.

	var_att_x=zeros(1,m+1);
	var_att_z=zeros(1,m+1);


%The main recursive loop.

	for k=0:m
	  GAIN=PMINUS*H'*inv(H*PMINUS*H'+R);
	  PPLUS=(eye(state_dim)-GAIN*H)*PMINUS;

	  %Add extra step here to assure symmetry of the P matrix.

	  PPLUS=(PPLUS+PPLUS')/2;


	  %Store desired output as elements of previously defined row
	  %vectors.  (Note that first element corresponds to k=0, next
	  %element to k=1,etc.  There will be a total of m+1 elements
	  %in the output vector.)

	  var_att_x(k+1)=PPLUS(3,3);
	  var_att_z(k+1)=PPLUS(9,9);

	  %Project ahead and symmetrize.

	  PMINUS=PHI11*PPLUS*PHI11'+ Q11;
	  PMINUS=(PMINUS+PMINUS')/2;
	end


	%Compute rms error vectors from variance vectors

	rms_x = sqrt(var_att_x);
	rms_z = sqrt(var_att_z);


%To view outputs while in workspace use appropriate MATLAB plot statements.
%Outputs can be saved with MATLAB save statements, if so desired.

%Scaling below is for MATLAB 4.  Scaling was avoided for MATLAB S3.5
%because of unresolved problems when involving semilog plots.
	semilogy(time,rms_x,time,rms_z);
	axis([0 1000 1e-4 1e-1]);

%Return graphics to auto scaling mode.
	clear axis;
