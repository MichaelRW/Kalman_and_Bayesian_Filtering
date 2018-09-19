function y=chisq16(x)
%CHISQ16   16 degree-of-freedom (DOF) chi-square density function.
%  See Problem 1.42 for general mathematical form.

y=(1/(256*5040)).*x.*x.*x.*x.*x.*x.*x.*exp(-.5.*x)
