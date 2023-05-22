function [orth,res1,res2,callcounts,dodotimes,eigtimes] ...
  = randomcluster(al,n,L,ep,numtrials)
%RANDOMCLUSTER   Helper function for SEC_4_1.

orth = nan(numtrials,1);
res1 = nan(numtrials,1);
res2 = nan(numtrials,1);
callcounts = nan(numtrials,1);
dodotimes = nan(numtrials,1);
eigtimes = nan(numtrials,1);
for trial = 1:numtrials
  x = -L*rand(n,1);
  y = -L*rand(n,1);
  lambda = cumsum(10.^x);
  mu = cumsum(10.^y);
  [U,~] = qr(randn(n)); U = U*diag(sign(randn(n,1)));
  E = randn(n); E = (E-E')/2;
  A = U*diag(lambda)*U';
  B = (U*expm(ep/2*E))*diag(mu)*(U*expm(ep/2*E))';
  dodotic = tic;
  [Uhat,Lambdahat,Mhat,callcount,eigtime] = dodo(A,B,al);
  dodotimes(trial) = toc(dodotic);
  eigtimes(trial) = eigtime;
  orth(trial) = norm(eye(n)-Uhat'*Uhat);
  normA = norm(A);
  normB = norm(B);
  res1(trial) = norm(Uhat'*A*Uhat-Lambdahat)/(normA+normB);
  res2(trial) = norm(Uhat'*B*Uhat-Mhat)/(normA+normB);
  callcounts(trial) = callcount;
end
orth = max(orth);
res1 = max(res1);
res2 = max(res2);
