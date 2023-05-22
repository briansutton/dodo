function [R,T] = grunbaum(n,r1,r2,r3,p)
% An example from "Toeplitz Matrices Commuting with Tridiagonal Matrices" by
% F. Albert Grunbaum.

if nargin<5
  p = acos(1/2*sqrt((r1^2-r1*r3)/(r2^2-r1*r3)));
end
r = [0;r1;r2;r3;arrayfun(@(k) r1*chebu(k-1,r2/r1*cos(p))./chebu(k-1,cos(p)),(4:n-1)')];
R = toeplitz(r);
N = n-1;
j = (1:n)';
a = 2*sin(j*p).*sin((j-N-1)*p);
b = 2*(r2-r1)/r1*cos(p).*cos(p*(2*j-N-2));
Dminus = spdiags([-ones(n,1) ones(n,1)],[0 -1],n,n);
Dplus = spdiags([-ones(n,1) ones(n,1)],[0 1],n,n);
T = diag(b)-Dminus*diag(a)*Dplus;
