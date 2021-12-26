%This file is for Problem 9.3, Part (b)

tau0=0.1;
for steps=0:10
    tau(steps+1)=tau0*2^steps;
    m=tau(steps+1)/tau0;
    sum=0;
    for i=1:N-2*m
        sum=sum+(xx(i+2*m)-2*xx(i+m)+xx(i))^2;
    end;
    sqrtav(steps+1)=sqrt(sum/2/(N-2*m)/tau(steps+1)^2);
end;
