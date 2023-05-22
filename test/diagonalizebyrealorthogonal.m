function [Q,La,stats] = diagonalizebyrealorthogonal(X)
%DIAGONALIZEBYREALORTHOGONAL   Diagonalize a complex symmetric normal matrix.
%
%   [Q,LA] = DIAGONALIZEBYREALORTHOGONAL(X) computes the eigenvalue
%   decomposition X = Q*LA*Q'. X must be complex symmetric and normal, and Q is
%   real orthogonal, not just complex unitary.
%
%   [...,STATS] = DIAGONALIZEBYREALORTHOGONAL(X) also calculates some accuracy
%   measurements for comparing this method with the built-in SCHUR.
%
%   Copyright 2023 Brian Sutton.

n = height(X);
offdiag = @(A) A-diag(diag(A));
stats(1) = norm(X-X.','fro');
stats(2) = norm(X'*X-X*X','fro');
[V,D] = schur(X);
stats(3) = norm(imag(V),'fro');
stats(4) = norm(eye(n)-V'*V,'fro');
stats(5) = norm(offdiag(V'*X*V),'fro');
A = real(X);
B = imag(X);
al = 0.01;
[Q,Lare,Laim] = dodo(A,B,al);
La = Lare+1i*Laim;
stats(6) = norm(imag(Q),'fro');
stats(7) = norm(eye(n)-Q'*Q,'fro');
stats(8) = norm(offdiag(Q'*X*Q),'fro');
