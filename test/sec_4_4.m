% Generates data as in Section 4.4 of the article "Simultaneous diagonalization
% of nearly commuting Hermitian matrices: do-one-then-do-the-other"

al = 0.01;
n = 7;
[j,k] = meshgrid(1:n,1:n);
A = cos((j+k)*pi/n);
A = A-diag(diag(A));
A = A+diag((2-n)/2*cos(2*(1:n)*pi/n));
B = sin((j+k)*pi/n);
B = B-diag(diag(B));
B = B+diag((2-n)/2*sin(2*(1:n)*pi/n));
[U,La,M] = dodo(A,B,al);
norm(eye(n)-U'*U,'fro')
norm(U'*A*U-La,'fro')
norm(U'*B*U-M,'fro')
