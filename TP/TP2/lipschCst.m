function L = lipschCst(op,opt,dim)
%% lIPSCHCST computes the Lipschitz constant of a given operator
%   L = lIPSCHCST(OP,OPT,IM) calculates the Lipschitz constant of operator
%   OP, with transpose OPT, in dimension DIM (DIM=1 for a signal, DIM=2 
%   for an image).
%
%
%   Example
%   ---------
% 
%   D   = @(x) opD(x);
%   Dt  = @(x) opDt(x);
%   lip = lipschCst(D,Dt,2)
% 
%   fprintf('The lipschitz constant %0g.\n',lip);
% 
%
%
%
%   H   = @(x) opH(x,'gaussian',5);
%   Ht  = @(x) opHt(x,'gaussian',5);
%   lip = lipschCst(H,Ht,1)
% 
%   fprintf('The lipschitz constant %0g.\n',lip);
%  


tol    = 1e-4;
lip(1) = 1;
rhon   = lip(1)+2*tol;

if      dim == 1, xn  = randn(64,1);
elseif  dim == 2, xn  = randn(64);
end
xnn = xn;

k = 1;
T = @(x) opt(op(x));

while abs(lip(k)-rhon)/lip(k) >= tol
   xn  = T(xnn);
   xnn = T(xn);

   rhon = lip(k);
   k    = k+1;
   lip(k) = sum(xnn(:).^2)/sum(xn(:).^2);
end

L = 1.01*(lip(k-1) + 1e-16);
%L = sqrt(L);


end




