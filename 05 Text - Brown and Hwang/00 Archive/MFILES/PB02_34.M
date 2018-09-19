clear
%File name is pb02_34.m and this file is for Problem 2.34.

%Part (a):  The program for generating the 256-sample realization is
%similar to that used in Problem 2.33* except for delta t.

%The rand statement used here is from MATLAB Version 3.5, and it
%results in a warning statement when using Version 4.0. This does
%not interfere with the solution.  If you are using Version 4.0,
%the warning statement can be eliminated by deleting the
%rand('normal') statement in line 24 of the code and replacing
%rand with randn everywhere.

varwk=1-exp(-2*1);
sigmawk=sqrt(varwk);
phi=exp(-1);

%Now set up 256-point time sample as a vector called gtime.

gtime=zeros(1,256);

%Next generate the 256-point sample realization.

rand('normal')
gtime(1)=rand;
for i=1:255
   gtime(i+1)=phi*gtime(i)+sigmawk*rand;
end

%Preview the sample realization for reasonableness.

plot(gtime)
title('Press ENTER to Continue')
pause

%Parts (b) and (c):  We want the periodgrams for the 64-point, 128-
%point, and 256-point samples of gtime.  We will first get the dft
%for these time samples and then form the periodograms from the dfts.
%dumxx's will be used as dummy variables.  (Note that the scaling
%here will not be the same as in the text.)

dum64=gtime(1:64);
dft64=fft(dum64);
abdft64=abs(dft64);
pgram64=abdft64.^2;

%Now plot the periodogram (first term is for zero frequency).

k=0:1:63;
plot(k,pgram64,'o')
title('Press ENTER to Continue')
pause

%Now use similar code and plot the 128-point periodogram.

dum128=gtime(1:128);
dft128=fft(dum128);
abdft128=abs(dft128);
pgram128=abdft128.^2;
kk=0:1:127;
plot(kk,pgram128,'o')
title('Press ENTER to Continue')
pause

%Now use similar code and plot the 256_point periodogram.

dft256=fft(gtime);
abdft256=abs(dft256);
pgram256=abdft256.^2;
kkk=0:1:255;
plot(kkk,pgram256,'o')
title('Press ENTER to Continue')
pause

%Part (d):  To do this part, begin with the periodogram for the
%256-point case.  The variable pgram256 is a 1 x 256 row vector,
%and we need to average successive blocks of 8 elements of this
%vector and then plot the results.  The plot will be truncated
%at the fold-over frequency.

avpgram=zeros(1,16);
for i=1:16
   avpgram(i)=(sum(pgram256((8*i-7):(8*i))))/8;
end

%For comparison purposes on the plot, the spacing between samples
%should be 8 times that of the 256-point plot. Therefore, set kkkk
%accordingly with the appropriate offset from zero, and truncate
%the right half of the plot.

kkkk=4:8:252;
truncpgm=[avpgram zeros(1,16)];
plot(kkkk,truncpgm,'o')
title('Press ENTER to end Problem 2.34')
