function y=gsq329(x)
%GSQ329   Special function to be integrated in Problem 3.29.
%  The inverse Laplace transform of G(s) is (1/w0)*sin(w0*t).
%  The scale factor for the input PSD is included in function y.

y=(10/(400*pi*pi)).*(sin(20*pi.*x)).*(sin(20*pi.*x))
