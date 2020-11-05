function Lpl = laguerre(p,l,X)
% LAGUERRE associated Laguerre function
%
% L = LAGUERRE(p,l,X) evaluate associated Laguerre function.
% p and l must be integer scalars greater than zero

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

x=X(:);

Lplt=zeros(numel(x),max([p+1,2]));

Lplt(:,1)=1;
Lplt(:,2)=-x+l+1;

for ii=2:p
    Lplt(:,ii+1)=1/ii*((2*ii+l-1-x).*Lplt(:,ii)-(ii+l-1).*Lplt(:,ii-1));
end

Lpl=reshape(Lplt(:,p+1),size(X));
