clear
%File name is ex11_1a.m and this is for Example 11.1.
%Initial PMINUS is modified in this version of the example.

%In this covariance analysis, each of the clock, selective avail-
%ability, and long-time-constant range error have 2-state models.
%The total system then has 6 states and the parameters are as given
%in the code below.

%The measurement sequence starts at k=0 and ends at k=m.

	m=18000;
	disp ('How many steps (default is 18000)');
	minput = input ('? ');
	if minput >= 0,
  	  m = minput;
	end;

%The noise processes are uncoupled, so the partioned 2 x 2 parts of
%phi and que will be computed separately, and then combined later.
%First compute phi and que's for the 3 clocks; sf and sg denote the
%input PSDs in the clock model.  (See Figure 11.5 of the text.)

dt=1
PHICLOCK=[1 dt;0 1];

sf1=.036;
sg1=.142;
QCLOCK1=[sf1*dt+sg1*(dt^3)/3 sg1*(dt^2)/2
         sg1*(dt^2)/2      sg1*dt];

sf2=.0143;
sg2=.0003;
QCLOCK2=[sf2*dt+sg2*(dt^3)/3 sg2*(dt^2)/2
         sg2*(dt^2)/2      sg2*dt];

sf3=.0018;
sg3=1.4e-10;
QCLOCK3=[sf3*dt+sg3*(dt^3)/3 sg3*(dt^2)/2
         sg3*(dt^2)/2      sg3*dt];

%Now compute phi and que for the SA model.  (SA model is the same as
%given in Problem 5.4*.  RMS value for SA is set at 23 m.)

w0=.012;
csq=.002585;
F=[0 1;-(w0^2) -(sqrt(2))*w0];
GWGT=[0 0;0 csq];
A=[-dt*F dt*GWGT;zeros(2) dt*F'];
B=expm(A);
PHISAT=B(3:4,3:4);
PHISA=PHISAT';
QSA=PHISA*B(1:2,3:4);

%Next, compute phi and que for long-time-constant 2nd-order process.
%Time constant is 30 min and its 2 repeated poles are at -1/1800  (on
%negative real axis).  RMS amplitude of the process is 5 meters.
%This process is referred to as the NONSA process.

w2=1/1800;
F=[0 1;-(w2^2) -2*w2];
GWGT=[0 0;0 100*(w2^3)];
A=[-dt*F dt*GWGT;zeros(2) dt*F'];
B=expm(A);
PHINONT=B(3:4,3:4);
PHINONSA=PHINONT';
QNONSA=PHINONSA*B(1:2,3:4);

%Now combine the partitioned parts of phi into a matrix called PHI.

PHI=[PHICLOCK zeros(2) zeros(2)
     zeros(2) PHINONSA zeros(2)
     zeros(2) zeros(2) PHISA];

%Now set up row matrices with the correct dimensions for your outputs.

	save_x=[zeros(1,m/60+1);
		zeros(1,m/60+1);
		zeros(1,m/60+1)];


%Start by specifying one of the three different sub-problems
%illustrated by this example.  (Each clock model is a sub-problem.)


      disp ('Running...');

      for mode=1:3

	if mode == 1,
	  disp ('  1. Compensated crystal');
	else if mode == 2,
	  disp ('  2. Ovenized crystal');
	  else if mode == 3,
	    disp ('  3. Rubidium standard');
	    end;
	  end;
	end;

	if mode == 1,
        Q=[QCLOCK1 zeros(2) zeros(2)
           zeros(2) QNONSA zeros(2)
           zeros(2) zeros(2) QSA];
        PMINUS = [1000000 0 0 0 0 0;
                  0 1000000 0 0 0 0;
                  0 0   25.0  0 0 0;
                  0 0   0  25/(1800*1800) 0 0;
                  0 0   0   0 529 0;
                  0 0   0   0 0 0.076176];
	else if mode == 2,
        Q=[QCLOCK2 zeros(2) zeros(2)
           zeros(2) QNONSA zeros(2)
           zeros(2) zeros(2) QSA];
        PMINUS = [1000000 0 0 0 0 0;
                  0 10000 0 0 0 0;
                  0 0   25.0  0 0 0;
                  0 0   0  25/(1800*1800) 0 0;
                  0 0   0   0 529 0;
                  0 0   0   0 0 0.076176];
	else if mode == 3,
        Q=[QCLOCK3 zeros(2) zeros(2)
           zeros(2) QNONSA zeros(2)
           zeros(2) zeros(2) QSA];
        PMINUS = [1000000 0 0 0 0 0;
                  0 1.0 0 0 0 0;
                  0 0   25.0  0 0 0;
                  0 0   0  25/(1800*1800) 0 0;
                  0 0   0   0 529 0;
                  0 0   0   0 0 0.076176];
     	end;
	end;
	end;
	state_dim = 6;


%Specify the numerical values of H and R for your problem.  (Caution:
%If either of these is a function of time, the assignment statement
%must be moved inside the recursive loop, i.e., the k loop.)

        H=[1 0 1 0 1 0; 0 1 0 1 0 1];
	R=[4.0 0; 0 0.002];


%Now set up row matrices with the correct dimensions for your outputs.

	var_x=zeros(1,m/60+1);
	

%Now do the main recursive loop.

	ksamp=0;
	for k=0:m
	  GAIN=PMINUS*H'*inv(H*PMINUS*H'+R);
	  PPLUS=(eye(state_dim)-GAIN*H)*PMINUS;

	%Add extra step here to assure symmetry of the P matrix.

	  PPLUS=(PPLUS+PPLUS')/2;


	%Now store desired output as elements of previously defined row vectors.
	%(Note that first element corresponds to k=0, next element to k=1,etc.
	%There will be a total of m+1 elements in the output vector.)

	  if rem(k,60)==0,
	    ksamp=ksamp+1;
	    var_x(ksamp)=PPLUS(1,1);
	  end;

	%Now project ahead and symmetrize.

	  PMINUS=PHI*PPLUS*PHI'+ Q;
	  PMINUS=(PMINUS+PMINUS')/2;
	end


	%Compute rms error vectors from variance vectors

	rms_x = sqrt(var_x);
	save_x(mode,:) = rms_x;

      end;

%To view outputs while in workspace use appropriate MATLAB plot statements.
%Outputs can be saved with MATLAB save statements.

	plot(save_x');
	xlabel('time (minutes)');
        title('Press ENTER to end Example 11.1')
