% Generates data as in Section 4.2 of the article "Simultaneous diagonalization
% of nearly commuting Hermitian matrices: do-one-then-do-the-other." Install
% the anymatrix library (https://github.com/mmikaitis/anymatrix) and add to
% your path before executing this code.

fmt = format('shorte');
[~,dims] = anymatrix('hadamard/complex_hadamard');
nk = [ 3 1; 4 1; 5 1; 6 1; 6 3; 7 1; 8 1; 9 1; 10 1; 11 2; 12 1; 13 1; 13 2 ];
complexhadamarddata = nan(height(nk),10);
for i = 1:height(nk)
  H = anymatrix('hadamard/complex_hadamard',nk(i,1),nk(i,2));
  [~,~,stats] = diagonalizebyrealorthogonal(H);
  complexhadamarddata(i,:) = [nk(i,1:2) stats];
end
complexhadamarddata(:,[1 2 5:10])
plusnoisedata = nan(height(nk),10);
for i = 1:height(nk)
  H = anymatrix('hadamard/complex_hadamard',nk(i,1),nk(i,2));
  G = randn(nk(i,1))+1i*randn(nk(i,1));
  G = (G+G.')/2;
  [~,~,stats] = diagonalizebyrealorthogonal(H+1e-6*G);
  plusnoisedata(i,:) = [nk(i,1:2) stats];
end
plusnoisedata(:,[1 2 5:10])
format(fmt);
