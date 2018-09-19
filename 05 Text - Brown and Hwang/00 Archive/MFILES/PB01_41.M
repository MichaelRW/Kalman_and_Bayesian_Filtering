clear
%File name is pb01_41.m and this file is for Problem 1.41.
                         
%We need to calculate the probability of k misses (i.e., errors)
%in n trials, where k is variable, i.e., k = 0, 1, 2, ... etc.
%The procedure is to calculate log(probability) first, and then get
%the probability as exp(log(probability)).
%First, specify tentative values for p, n, and k.

p=.001
n=5000
k=5

%Next, set up a vector for the probabilities of exactly 0, 1, 2, ... k
%misses, and calculate the probability of exactly 0 misses separately.
%Note that the first element of vector prob corresponds to 0 misses.

prob=zeros(1,k+1)
dum1=n*log(1-p)
prob(1)=exp(dum1)

%Next, work out the remaining k elements of prob by breaking up log(p)
%into 4 additive parts.  (You will need to do some paper-and-pencil
%algebraic reduction before doing this step.)

a=0
b=0
for i=1:k
   a=a+log(n+1-i)
   b=b+log(i)
   c=i*log(p/(1-p))
   d=n*log(1-p)
   dum2=a-b+c+d
   prob(i+1)=exp(dum2)
end
cumprob=sum(prob)

%Now, using your editor, iterate on k until you find a value which
%will yield a cumulative probability (i.e., cumprob) of about .98.
%Ans:  k=10
