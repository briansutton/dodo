% Generates data as in Section 4.1 of the article "Simultaneous diagonalization
% of nearly commuting Hermitian matrices: do-one-then-do-the-other"

al = 0.01;
numtrials = 100;
n = [10 10 100 100 1000 1000 10 10 100 100 1000 1000]';
L = [5 20 5 20 5 20 5 20 5 20 5 20 ]';
ga = [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3]';
orth = nan(size(n));
res1 = nan(size(n));
res2 = nan(size(n));
callcounts = nan(numtrials,12);
dodotimes = nan(numtrials,12);
eigtimes = nan(numtrials,12);
for k = 1:length(n)
  k
  [orth(k),res1(k),res2(k),callcounts(:,k),dodotimes(:,k),eigtimes(:,k)] = ...
    randomcluster(al,n(k),L(k),ga(k),numtrials);
end
fmt = format('shorte');
[ n L ga res1 res2 orth res1./ga res2./ga ]
format(fmt);
figure;
timedata = reshape(dodotimes(:,[5 6 11 12]),4*numtrials,1);
histogram(timedata,linspace(0,round(2*max(timedata),1,'significant'),40));
hold on;
medeigtime = eigtimes(:,[5 6 11 12]); medeigtime = median(medeigtime(:));
xline(medeigtime,'k');
xlabel('running time (s)');
ylabel('frequency');
yticks([]);
figure;
calldata = reshape(callcounts(:,[5 6 11 12]),4*numtrials,1);
histogram(calldata,linspace(0,round(1.33*max(calldata),1,'significant'),40));
xlabel('invocations of DODO');
ylabel('frequency');
yticks([]);
