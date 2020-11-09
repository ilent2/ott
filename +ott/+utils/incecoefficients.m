function [A_n,B_n,eigA,eigB] = incecoefficients(p,xi)
%INCECOEFFICIENTS calculates fourier coefficients for Ince polynomials
%
% [A,B,eigA,eigB] = INCECOEFFICIENTS(p,xi) calculate the two sets of
% fourier coefficients of Ince polynomials for a given order, p,
% and ellipticity, xi, in eigen-matrix form.
%
% Notes:
% A_n is the eigen-matrix corresponding to the form:
% \Sum_i A(n,i)*cos((2*(ii-1)+delta)*z) where delta=1 for p odd and 0 for
% even. 1<=i<=p.
%
% B_n is the eigen-matrix corresponding to the form:
% \Sum_i B(n,i)*sin((2*(ii-1)+delta)*z) where delta=1 for p odd and 0 for
% even. 1<=i<=p.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

delta=ceil(p/2)+1~=p/2+1;
W1=zeros(ceil(p/2)+1,ceil(p/2)+2);
W1(1,1:2)=[0,xi*(p/2+1)];
W1(2,1:3)=[xi*p,4*(1-delta)+delta*(1+xi/2*(p+1)),xi/2*(p+4-delta)];
for ii=2:size(W1,2)-2
    W1(ii+1,ii:ii+2)=[xi/2*(p-2*ii+2+delta), ...
        (2*ii-delta)^2,xi/2*(p+2*ii+2-delta)];
end

W1(:,end)=[];

[A_n,eigsA]=eig(W1(delta+1:end,delta+1:end));
ev=diag(eigsA);
ev(ev==0&any(A_n==1,1)')=-inf; % this checks for zero eigenvalue problems
[valz,indx]=sort(ev);             % <\ not necessary if done properly

valz(isinf(valz))=0;
eigA=valz;

permMat=zeros(size(A_n));
permMat(((1:size(A_n,1))-1).'*size(A_n,1)+indx)=1; % we need to permute
A_n=A_n*permMat;          % <\ the eigenvectors from lowest to highest

b=sign(A_n(1+sum(valz(1)==0),:));
b(b==0)=1;
A_n=repmat(b,[size(A_n,1),1]).*A_n;
A_n=A_n'; %we like row vectors

if delta
    W1(2,2)=1-xi/2*(p+1);
end

[B_nt,eigsBt]=eig(W1(2:end,2:end));

B_n=eye(size(A_n));
eigsB=zeros(size(A_n));

B_n(2-delta:end,2-delta:end)=B_nt;
eigsB(2-delta:end,2-delta:end)=eigsBt;

ev=diag(eigsB);
ev(ev==0&any(B_n==1,1)')=-inf; % this checks for zero eigenvalue problems.
[valz,indx]=sort(ev);             % \> not necessary if done properly

valz(isinf(valz))=0;
eigB=valz;

permMat=zeros(size(B_n));
permMat(((1:size(B_n,1))-1).'*size(B_n,1)+indx)=1; % we need to permute
B_n=B_n*permMat;            % \> the eigenvectors from lowest to highest

b=sign(B_n(1+sum(valz(1)==0&p~=1),:));
b(b==0)=1;
B_n=repmat(b,[size(B_n,1),1]).*B_n;

B_n=B_n'; %we like row vectors


