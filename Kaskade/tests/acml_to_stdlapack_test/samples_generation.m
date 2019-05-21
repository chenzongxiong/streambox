% Generation of test examples for the linalg/acml_to_lapack.hh c++
% headerfile

format long
A=[1 2 3; 4 5 6; 7 8 10; 10 11 14; 13 16 25];
B=[9 6 3; 8 5 2; 7 4 1];
C=[-6 -5 -4; -9 -8 -7; -3 -2 -1; 1 2 3; 4 5 6];
alpha=3;
beta=4;
fprintf(1,'dgemm:');
C=alpha*A*B'+beta*C
x=[1 4 9 16 25]'
y=[9 4 1]';
fprintf(1,'dgemv:'); 
y=alpha*A'*x + beta*y
% print A
AQR=A(1:3,1:3);
x=[1 4 9]';
fprintf(1,'dgesv:'); 
z=linsolve(AQR,x)
AQR=A;
x=[1 4 9 16 25]';
fprintf(1,'dgels,dgelsy:');
z=lsqlin(A, x)
AQR=A(1:3,1:3);
fprintf(1,'dgetrf&dgetri,sgetrf&sgetri:');
Ainv=inv(AQR)
AQR=[ 1 2 3; 2 3 5; 3 5 7];
fprintf(1,'dsyevd:');
[V,D]=eig(AQR)
ASVD=A;
fprintf(1,'dgesvd:');
[U,S,V] = svd(ASVD)
V=V'
AQR=[ 1 2 0; 2 3 5; 0 5 7];
fprintf(1,'dstev:');
[V,D]=eig(AQR)
V=V'

