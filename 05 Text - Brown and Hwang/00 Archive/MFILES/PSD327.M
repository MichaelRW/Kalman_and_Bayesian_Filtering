function y=psd327(x)
%PSD327   Special function to be integrated in Problem 3.27.

y=((1./pi).*((sin(.5.*x))./(.5.*x)).^2)./(1+x.*x)
