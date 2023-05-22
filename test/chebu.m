function y = chebu(n,x)
%CHEBU   ChebyshevU polynomial.

if n == 0, y = 1; return; end
a = 1;
b = 2*x;
k = 1;
for k = 2:n
  c = 2*x*b-a;
  a = b;
  b = c;
end
y = b;
