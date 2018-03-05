function [modeweights,LGlookups,IGlookups]=genLG2IG(order_paraxial,xi)
% genLG2IG.m --- LG->IG conversion matrix for elipticity xi. A 
%		parameter of xi=0 gives the non-vortex LG modes.
%		a parameter of xi=1e100 will give pretty HG modes.
%
% Usage:
%
% [modewieghts,LGlookups,IGlookups] = genLG2IG(paraxial_order,xi)
%
% PACKAGE INFO

%first create the upper block... these are the fourier coefficients...
[A_n,B_n]=incecoefficients(order_paraxial,xi);

%prepare the index matrices for this upper block.
p=[floor(order_paraxial/2):-1:0];
l=order_paraxial-2*p;

P=repmat(p,[length(p),1]);
L=repmat(l,[length(l),1]);

Nb=(sqrt(factorial(P+L).*factorial(P)).*(-1).^(P+L+(order_paraxial+L.')/2));
Na=(1+1*(L==0)).*(sqrt(factorial(P+L).*factorial(P)).*(-1).^(P+L+(order_paraxial+L.')/2));

NA_n=flipud(Na.*A_n);
NB_n=flipud(Nb.*B_n);

%calculate the other blocks:
if rem(order_paraxial,2)==1
    bigA=[fliplr(NA_n),NA_n];
    bigB=[fliplr(NB_n),-NB_n];

    modeweights=reshape([bigA.';bigB.'],[order_paraxial+1,order_paraxial+1]).';
    
else
    bigA=[fliplr(NA_n),NA_n(:,2:end)];
    bigB=[fliplr(NB_n),-NB_n(:,2:end)];
    
    modeweights=reshape([bigA.';bigB.'],[order_paraxial+1,order_paraxial+2]).';
    modeweights(end,:)=[];
end

for ii=1:size(modeweights,1)
    modeweights(ii,:)=modeweights(ii,:)/sqrt(sum(abs(modeweights(ii,:)).^2));
end

%create imaginary matrix:
imat=repmat((-1i).^[0:order_paraxial].',[1,order_paraxial+1]);

modeweights=modeweights.*imat;

if nargout>1
	[LGlookups,IGlookups]=meshgrid([0:order_paraxial],[0:order_paraxial]);
end
