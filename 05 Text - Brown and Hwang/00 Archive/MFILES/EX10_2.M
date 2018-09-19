clear
%File name is ex10_2.m and this is for Example 10.2 (INS/DME)

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

%The measurement sequence starts at k=0 and ends at k=m.

	m=200;
	disp ('How many steps (default is 200)');
	minput = input ('? ');
	if minput >= 0,
  	  m = minput;
	end;
	

%Now set up row matrices with the correct dimensions for your outputs.

	save_x=[zeros(1,m+1);
		zeros(1,m+1);
		zeros(1,m+1)];
	save_y=[zeros(1,m+1);
		zeros(1,m+1);
		zeros(1,m+1)];
	time=0:1:m;


%Start by specifying one of the three different sub-problems
%illustrated by this example.

      disp ('Running...');

      for mode=1:3

	if mode == 1,
	  disp ('  1. INS-only');
	else if mode == 2,
	  disp ('  2. INS/DME');
	  else if mode == 3,
	    disp ('  3. DME-only');
	    end;
	  end;
	end;



%Specify the initial P matrix as PMINUS.

	if mode <= 2,
	  insphiq;
          PMINUS= [10 0  0  0  0  0  0  0  0;
                   0  0.000001 0  0  0  0  0  0  0;
                   0  0  0.000001 0  0  0  0  0  0;
                   0  0  0  10 0  0  0  0  0;
                   0  0  0  0  0.000001 0  0  0  0;
                   0  0  0  0  0  0.000001 0  0  0;
                   0  0  0  0  0  0  10 0  0;
                   0  0  0  0  0  0  0  0.000001 0;
                   0  0  0  0  0  0  0  0  0.000001];
     	  state_dim = 9;
	else if mode == 3,
       	    PHI = [1 0; 0 1];
       	    Q = [400 0; 0 400];
       	    PMINUS = [10 0; 0 10];
       	    state_dim = 2;
     	  end; 
	end;


%Specify the numerical values of H and R for your problem.  (Caution:
%If either of these is a function of time, the assignment statement
%must be moved inside the recursive loop, i.e., the k loop.)

	if mode <= 2,
     	  H=[zeros(1,9); zeros(1,9)];
	else
     	  H=[0 0; 0 0];
	end;

	R=[225 0; 0 225];


%Now set up row matrices with the correct dimensions for your outputs.

	var_x=zeros(1,m+1);
	var_y=zeros(1,m+1);

     
%Now do the main recursive loop.

	for k=0:m
	  if mode >= 2,
	    drncos_x=10000;
	    drncos_y=-10000+100*k;
	    magnitude=sqrt (drncos_x^2+drncos_y^2);
	  end;
	
	  if mode == 2,
 	    H(1,1)=drncos_x/magnitude;
	    H(1,4)=drncos_y/magnitude;
	    H(2,1)=-drncos_x/magnitude;
	    H(2,4)=drncos_y/magnitude;
          else if mode == 3,
	      H(1,1)=drncos_x/magnitude;
	      H(1,2)=drncos_y/magnitude;
	      H(2,1)=-drncos_x/magnitude;
	      H(2,2)=drncos_y/magnitude;
	    end;
	  end;
	  GAIN=PMINUS*H'*inv(H*PMINUS*H'+R);
	  PPLUS=(eye(state_dim)-GAIN*H)*PMINUS;	
 	
	%Add extra step here to assure symmetry of the P matrix.

	  PPLUS=(PPLUS+PPLUS')/2;


	%Now store desired output as elements of previously defined row vectors.
	%(Note that first element corresponds to k=0, next element to k=1,etc.
	%There will be a total of m+1 elements in the output vector.)

	  var_x(k+1)=PPLUS(1,1);
	  if state_dim == 2,
	    var_y(k+1)=PPLUS(2,2);
	  else
	    var_y(k+1)=PPLUS(4,4);
	  end;


	%Now project ahead and symmetrize.

	  PMINUS=PHI*PPLUS*PHI'+ Q;
	  PMINUS=(PMINUS+PMINUS')/2;
	end


	%Compute rms error vectors from variance vectors

	rms_x = sqrt(var_x);
	rms_y = sqrt(var_y);
	save_x(mode,:) = rms_x;
	save_y(mode,:) = rms_y;
      end;

%To view outputs while in workspace use appropriate MATLAB plot statements.
%Outputs can be saved with MATLAB save statements, if so desired.
      
	semilogy(time,save_x(1,:),time,save_x(2,:),time,save_x(3,:));
        title('Press ENTER to Continue')
	pause
	semilogy(time,save_y(1,:),time,save_y(2,:),time,save_y(3,:));
	title('Press ENTER to end Example 10.2')
     
