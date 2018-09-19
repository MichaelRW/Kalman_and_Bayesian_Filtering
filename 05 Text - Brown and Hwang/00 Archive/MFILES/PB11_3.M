clear
% File name is PB11_3.M
% The purpose of this m-file is to solve for the state transition matrix 
% and the Q matrix in Problem 11.3 using the symbolic functions of Matlab
% Version 4.

% Declaring the F and G matrices; sigma, beta.

F=sym('[0 1 0;0 0 1;0 0 -beta]')
G=sym('[0; 0; (2*sigma^2*beta)^0.5]')
state_dim=3;


% The symbolic toolbox in Matlab 4 does not yet allow the 
% inverse Laplace transform function, invlaplace, to evaluate
% matrices.  Instead, the invlaplace function must evaluate
% each component of the matrix individually to obtain the
% state transition matrix STM.

sI=symmul('s',eye(state_dim))
sI_F=symsub(sI,F)
inv_sI_F=inverse(sI_F)
STM=zeros(state_dim,state_dim)
for i=1:state_dim
  for j=1:state_dim
    STM=sym(STM,i,j,invlaplace(sym(inv_sI_F,i,j)))
  end
end;

% The double integral of Eq. 5.3.6 must be evaluated to obtain Q.
% However, inner integral may be implicitly evaluated using the 
% sifting property of the Dirac delta function.  We only need to
% deal with the outer integral where the symbolic toolbox allows 
% matrix integration.

GWGt=symmul(G,transpose(G))
INT23=symmul(GWGt,transpose(STM))
INT123=symmul(STM,INT23)
Q=int(INT123,'t',0,'dt')
simplify(Q)
