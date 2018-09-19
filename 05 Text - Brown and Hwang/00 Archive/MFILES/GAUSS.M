function y=gauss(x,mn,sigma)
%GAUSS   General normal density function.
%  The first argument is the usual independent variable for
%  the function, just as with any of the single-variable
%  MATLAB built-in functions.  The second argument is the mean
%  (a scalar), and the third argument is the standard devi-
%  ation (a scalar).  When calling gauss, the arguments must
%  appear in the proper order (separated by commas) within 
%  the parentheses.

expcoeff=1/(2*sigma^2);
coeff=1/(sigma*sqrt(2*pi));
y=coeff.*exp(-expcoeff.*(x-mn).*(x-mn));
