function [output]=hermite(n,X);
%HERMITE calculate hermite polynomial.
%
% H = HERMITE(n,X) calculates hermite polynomial of order n and argument X.
% n must be integer scalars greater than zero

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

x=X(:);

h=zeros(numel(x),max([n+1,2]));

h(:,1)=1;
h(:,2)=2*x;

for ii=2:n
    h(:,ii+1)=2*x.*h(:,ii)-2*(ii-1)*h(:,ii-1);
end

output=reshape(h(:,n+1),size(X));
