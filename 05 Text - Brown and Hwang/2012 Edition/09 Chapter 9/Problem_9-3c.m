%This file is for Problem 9.3, Part (c)

h0=Sf*2;
h_2=Sg/2/(pi^2);
 
i=0;
for k=-1:0.02:2
    i=i+1;
    taux(i)=10^k;
    avar1_approx=h0/2/taux(i);
    avar3_approx=(2*pi)^2*h_2*taux(i)/6;
    adevasymp_approx(i)=max([avar1_approx,avar3_approx])^0.5;
end;
 
loglog(tau,sqrtav,'*-',taux,adevasymp_approx,'--');
axis([0.1 100 2E-12 2E-10]);
grid;
xlabel('Averaging Time (s)');
ylabel('Allan deviation');
legend('Empirical Allan Deviation','Theoretical Asymptotes');
