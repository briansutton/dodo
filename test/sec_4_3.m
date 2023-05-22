% Generates data as in Section 4.3 of the article "Simultaneous diagonalization
% of nearly commuting Hermitian matrices: do-one-then-do-the-other"

n = 100;
al = 0.1;
be = 0.2;
[j,i] = meshgrid(1:n,1:n);
R = sin(al*(i-j))./sin(be*(i-j));
R(i==j) = al/be;
[~,T] = grunbaum(n,R(1,2),R(1,3),R(1,4),be);
[VR,DR] = eig(R);
[~,index] = sort(diag(DR));
VR = VR(:,index); DR = DR(index,index);
[VT,DT] = eig(T);
[~,index] = sort(diag(DT));
VT = VT(:,index); DT = DT(index,index);
[U,La,M] = dodo(R,T);
offdiag = @(A) A-diag(diag(A));
fmt = format('shorte');
[norm(offdiag(VR'*R*VR),'fro') norm(offdiag(VR'*T*VR),'fro');
 norm(offdiag(VT'*R*VT),'fro') norm(offdiag(VT'*T*VT),'fro')]
[norm(U'*R*U-La,'fro') norm(U'*T*U-M,'fro')]
format(fmt);
