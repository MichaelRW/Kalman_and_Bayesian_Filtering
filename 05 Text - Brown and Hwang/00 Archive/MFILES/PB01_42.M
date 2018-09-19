clear
%File name is pb01_42.m and this file is for Problem 1.42.

%You will need the special functions for n = 1, 2, 4, and 16 degrees-
%of-freedom (DOF) chi-square densities for plotting.  These are
%provided in files chisq1.m, chisq2.m, chisq4.m, chisq16.m.

%You will also need the general normal density function with mean
%mn and standard deviation sigma, and this is provided as a MATLAB
%function function gauss(x,mn,sigma).  Use of this function is simi-
%lar to the usual built-in MATLAB functions such as sin(x), cos(x),
%etc, except that the 2nd and 3rd arguments (which are scalars)
%must also be specified when calling gauss(x,mn,sigma).  The 1st
%argument x is, in general, a vector and it is treated the same
%here as in any of the built-in single variable functions.

%This problem is a plotting exercise.  Try these statements 
%for part (a):

     x=.25:.25:10;
     plot(x,chisq1(x),x,chisq2(x),x,chisq4(x))
     title('Press ENTER to Continue')
     pause

%A scale change is needed for part (b).  Try:

     x=0:.5:40;
     mn=16
     sigma=sqrt(32)
     plot(x,gauss(x,mn,sigma),x,chisq16(x))
     title('Press ENTER to end Problem 1.42')
