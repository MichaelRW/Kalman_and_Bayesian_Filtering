%This file is for Problem 7.7

function [x,w] = resample(x,w,N)
 
u = rand(N,1);
wc = cumsum(w);
wc = wc/wc(N);
[dum,indl] = sort([u;wc]);
ind2 = find(indl <= N);
ind = ind2 - (0:N-1)';
 
x=x(ind,:);
w=ones(1,N)./N;
