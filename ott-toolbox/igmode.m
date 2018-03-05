function [IGo,IGe]=igmode(p,m,xi,x,y,z,theta);
% igmode.m --- Generates an Ince Gauss mode of order p, m and xi. xi is the
% ellipticity at z=0.
%
% Usage:
% [IGo,IGe] = igmode(p,m,xi,x,y);
%
% Notes:
% we use x,y input because the coordinate transformation has a hard cap as
% xi->1/eps.
%
% Given that we also do transformations if m is in the form [m,parity].
% it will only output EITHER the odd or even form parity=1 for the even
% form.
%
% PACKAGE INFO

if nargin<7
    theta=atan(1/pi)*180/pi;
end

if nargin<6
    z=0;
end

%we need to do this becasue of a numerical difficulty...
if xi>1e9
    xi=1e9;
end

w = 1.; %Beam waist in normalized units.

k=2*pi;
w0=w/pi/tan(abs(theta/180*pi));
zr=pi*w0^2;
w_z=w0*sqrt(1+(z/zr).^2);
R=1./z*(z.^2+zr.^2);
R(z==0)=inf;
psi=atan2(z,zr);

%the u and v we have is for z... resultingly we need to rescale x and y.
%and recalculate u and v. This is inelegant... but it works. Maybe we
%should have x, y input?

% xi is the parameter when w_z=w0; f_0=xi/2*w0^2

r=sqrt(x.^2+y.^2);

extFac=exp(-1i*(k*z+r.^2/2./R));

X=x./w_z;
Y=y./w_z;

[u,v]=xy2uv(X,Y,sqrt(xi*w0^2/2)); %sqrt(xi*w0^2/2)

r=r./w_z;

if p>=2
    delta=ceil(p/2)+1~=p/2+1;
    
    [A_n,B_n]=incecoefficients(p,xi);
    
    %% normalisation from here
    [M,P]=meshgrid([1:floor(p/2)+1],ones(size(A_n(:,1))));
    
    normmat=ones(size(A_n));
    
    if ~delta
        normmat(:,1)=2;
        normmat=normmat;
    end
    
    anorm=sqrt(1/2*abs(sum(normmat.*A_n.^2.*gamma((p-(2*(M-1)+delta))/2+1).*gamma((p+2*(M-1)+delta)/2+1),2)));
    bnorm=sqrt(1/2*abs(sum(B_n.^2.*gamma((p-(2*(M-1)+delta))/2+1).*gamma((p+2*(M-1)+delta)/2+1),2)));
    
    A_n=A_n./repmat(anorm,[1,size(A_n,2)]);
    B_n=B_n./repmat(bnorm,[1,size(A_n,2)]);
    
    %% Set up initial values
    C2=zeros(size(r));
    C1=C2;
    S1=C1;
    S2=C1;
    
    % these will be the normalisation constants.
    c0=0;
    cp=0;
    s0=0;
    sp=0;
    
    A_0=A_n(floor(m(1)/2)+1,1);
    B_0=B_n(floor(m(1)/2)+1,2-delta); %special case.
    
    for ii=1:size(A_n,2)
        C1=C1+A_n(floor(m(1)/2)+1,ii)*cosh((2*(ii-1)+delta)*u);
        C2=C2+A_n(floor(m(1)/2)+1,ii)*cos((2*(ii-1)+delta)*v);
        c0=c0+A_n(floor(m(1)/2)+1,ii);
        if delta
            cp=cp-A_n(floor(m(1)/2)+1,ii)*(cos((2*(ii-1)+delta)*pi/2)-cos((2*(ii-1)+delta)*(pi/2+1e-10)))/1e-10;
        else
            cp=cp+A_n(floor(m(1)/2)+1,ii)*cos((2*(ii-1)+delta)*pi/2);
        end
    end
    
    for ii=1:size(B_n,2)
        S1=S1+B_n(floor(m(1)/2)+1,ii)*sinh((2*(ii-1)+delta)*u);
        S2=S2+B_n(floor(m(1)/2)+1,ii)*sin((2*(ii-1)+delta)*v);
        
        s0=s0+(B_n(floor(m(1)/2)+1,ii)*sin((2*(ii-1)+delta)*1e-10))/1e-10;
        if delta
            sp=sp+B_n(floor(m(1)/2)+1,ii)*sin((2*(ii-1)+delta)*pi/2);
        else
            sp=sp-B_n(floor(m(1)/2)+1,ii)*(sin((2*(ii-1)+delta)*pi/2)-sin((2*(ii-1)+delta)*(pi/2+1e-10)))/1e-10;
        end
    end
    
    if m(1)==0
        B_0=1;
        sp=1;
        s0=1;
    end
    
    S1S2=S1.*S2/sp/s0*B_0*gamma((p+2-delta)/2+1)*sqrt(xi)*sqrt(delta+(1-delta)*xi); %
    C1C2=C1.*C2/cp/c0*A_0*gamma((p+delta)/2+1)*sqrt(xi)/sqrt(delta+(1-delta)*xi); %
else
    S1S2=1;
    C1C2=1;
    if p==1
        S1S2=sinh(u).*sin(v)*sqrt(xi*2);
        C1C2=cosh(u).*cos(v)*sqrt(xi*2);
    end
    
end

%calculate normalisations
IGe = 1i^(p+delta)*sqrt(2/pi)*exp(1i*(p + 1)*psi).*C1C2.*exp( -r.^2).*extFac./w_z;
IGo = 1i^(p+delta)*sqrt(2/pi)*exp(1i*(p + 1)*psi).*S1S2.*exp( -r.^2).*extFac./w_z;

if length(m)>1
	if m(2)==1
		IGo=IGe;
	end
end
