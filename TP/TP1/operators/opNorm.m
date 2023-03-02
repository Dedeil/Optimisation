function N = opNorm(op,opt,dim)
%% OPNORM computes the Lipschitz constant of a given operator
%   L = OPNORM(OP,OPT,IM) calculates the Lipschitz constant of operator
%   OP, with transpose OPT, in dimension DIM (DIM=1 for a signal, DIM=2 
%   for an image).
%
%
%   Example
%   ---------
% 
%   D  = @(x) opD(x);
%   Dt = @(x) opDt(x);
%   n  = opNorm(D,Dt,2)
% 
%   fprintf('The operator norm is: %0g.\n',n);
% 
%
%
%
%   H  = @(x) opH(x,'gaussian',5);
%   Ht = @(x) opHt(x,'gaussian',5);
%   n  = opNorm(H,Ht,1)
% 
%   fprintf('The operator norm is: %0g.\n',n);
%  


tol     = 1e-4;
norm(1) = 1;
rhon    = norm(1)+2*tol;

if      dim == 1, xn  = randn(64,1);
elseif  dim == 2, xn  = randn(64);
end
xnn = xn;

k = 1;
T = @(x) opt(op(x));

while abs(norm(k)-rhon)/norm(k) >= tol
   xn  = T(xnn);
   xnn = T(xn);

   rhon = norm(k);
   k    = k+1;
   norm(k) = sum(xnn(:).^2)/sum(xn(:).^2);
end

N = 1.01*(norm(k-1) + 1e-16);
%L = sqrt(L);


end




