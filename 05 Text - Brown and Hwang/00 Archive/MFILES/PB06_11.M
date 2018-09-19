clear
%File name is pb06_11.m and this is for Problem 6.11.

%The steady-state variables H, PHI, and gain were saved from Problem 6.10
%in a mat-file called prob610.mat.  These variables are needed here.
%(If you have not run pb06_10.m, do so before running this problem.)

load prob610

%Form the matrix for which we want the eigenvalues.

A=PHI-gain*H*PHI

%Now look at the eigenvalues of A.

b=eig(A)

%The error covariance in Prob. 6.10 approaches steady-state in about
%100 steps.  This is consistent (qualitatively) with the eigenvalues.
%The controlling eigenvalue is .9631, and this corresponds to a time
%constant of 26.6 steps.  Four time constants would be about 100 steps.
