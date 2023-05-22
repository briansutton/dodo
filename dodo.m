function [U,La,M,numcalls,eigtime] = dodo(A,B,al)
%DODO   Simultaneously diagonalize commuting Hermitian matrices.
%
%   [U,LA,M] = DODO(A,B) simultaneously diagonalizes commuting Hermitian
%   matrices A and B. It computes a unitary matrix U and diagonal matrices LA
%   and M for which U'*A*U = LA and U'*B*U = M. If A and B do not exactly
%   commute, then the residuals LA-U'*A*U and M-U'*B*U are nonzero.
%
%   [...] = DODO(A,B,AL) uses tolerance AL (alpha), 0 < AL <= 1, for
%   determining eigenvalue clusters. Smaller AL trades faster execution for
%   lower accuracy. Default: AL = 0.01.
%
%   [U,LA,M,NUMCALLS,EIGTIME] = DODO(...) also returns the number of
%   invocations of this recursive method and the execution time for the largest
%   eigenvalue-eigenvector computation.
%
%   Copyright 2023 Brian Sutton

if nargin<3, al = 0.01; end
n = size(A,1);
sAest = 1/sqrt(n)*norm(A-trace(A)/n*eye(n),'fro');
sBest = 1/sqrt(n)*norm(B-trace(B)/n*eye(n),'fro');
if sAest>=sBest
  [U,La,M,numcalls,eigtime] = dodocore(A,B,al,n);
else
  [U,M,La,numcalls,eigtime] = dodocore(B,A,al,n);
end

end



function [U,La,M,numcalls,eigtime] = dodocore(A,B,al,n)

numcalls = 1;
eigtic = tic; [U,La] = eig(hermpart(A)); eigtime = toc(eigtic);
la = diag(La);
[la,I] = sort(la); La = La(I,I); U = U(:,I);
sA = la(end)-la(1);
clusters = [ 1; 1+find(diff(la)>=(al/n)*sA) ];
A = diag(la);
B = U'*B*U;
for k = 1:length(clusters)
  i = clusters(k);
  if k==length(clusters), j = n; else, j = clusters(k+1)-1; end
  if j>i
    [Uk,Lak,Mk,numcallsk] = dodo(A(i:j,i:j),B(i:j,i:j),al);
    U(:,i:j) = U(:,i:j)*Uk;
    A(i:j,i:j) = Lak;
    B(i:j,i:j) = Mk;
    numcalls = numcalls+numcallsk;
  end
end
La = diag(diag(A));
M = diag(diag(B));

end



function A = hermpart(A)

n = size(A,1);
for i = 1:n
  for j = 1:i-1
    A(i,j) = (A(i,j)+conj(A(j,i)))/2;
    A(j,i) = A(i,j);
  end
  A(i,i) = real(A(i,i));
end

end
