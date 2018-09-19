function y=chisq1(x)
%CHISQ1   One degree-of-freedom (DOF) chi-square density function.
%  See Problem 1.42 for general mathematical form.

y=(1/sqrt(2*pi)).*((x).^(-.5)).*exp(-.5.*x)
