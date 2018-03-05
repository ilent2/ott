function [nn,mm,a,b] = bsc_lgmode_farfield( anglez,lg_mode,varargin )
% bsc_lgmode_farfield.m
% Uses a line integral method to find the
% spherical harmonic expansion of a Laguerre-Gaussian laser beam.
%
% Usage:
% [n,m,a,b] = bsc_lgmode_farfield( angle,lg_mode );
%
% calculates the expansion coefficients for a 'lg_mode' ([p,l]) x-polarised
% beam converging with a cone of 'angle' (in degrees)
%
% OR
%
% [n,m,a,b] = bsc_lgmode_farfield( angle,lg_mode,(optional))
%
% (optional) inputs are not limited to a single entry or order and include:
%
% truncation_angle --- angle of hard edge aperture (degrees) default: 90
% polarisation --- either a Jones vetor {[E_x,Ey] | 'azimuthal' |
% 'radial'}. default: [1,0]
% w0_scaling --- {'fixed' | {'free'} }. 'fixed' uses the same w0 for all modes,
% for mode matching this is essential. 'free' moves the waist for optimal OAM.
% NYI!!!: mapping --- {'direct' | {'projection'}}. 'direct' uses non-paraxial
% mapping. 'projection' uses the symmetry maintinging
% paraxial-to-non-paraxial convention.
%
% PACKAGE INFO

truncation_angle=90;
polarisation=[1,0];
mapping='projection';
w0_scaling='free';
paraxial_order=2*lg_mode(1)+abs(lg_mode(2));

tol=5e-7; %tol to truncate nmax;
itlim=200; %we have to stop somewhere!

debug=0;

%setting up
if nargin>2
    
    for ii=1:length(varargin)
        
        switch class(varargin{ii})
            case {'double'}
                if length(varargin{ii})>1
                    polarisation=varargin{ii};
                else
                    truncation_angle=varargin{ii};
                end
            case {'char'}
                switch lower(varargin{ii})
                    case 'azimuthal'
                        polarisation='az';
                    case 'radial'
                        polarisation='ra';
                    case 'direct'
                        mapping = 'direct';
                    case 'fixed'
                        w0_scaling='fixed';
                end
        end
    end
end

w = 1.; %Beam waist in normalized units.

if and(paraxial_order ~= 0 , strcmp(w0_scaling,'free'))
    invL=1./abs(paraxial_order );
    z = exp(-(abs(paraxial_order )+2.)*invL);
    w=-(1.+2*sqrt(invL)+invL); %This is a really good starting guess. It converges within 3 iterations for l=1:10000+
    
    w0=-w;
    
    while (abs(w-w0)>0.00001)
        w0=w;
        expw = exp(w);
        
        w=w0-(w0*expw+z)/(expw+w0*expw); %Newton's rule... Usually this kind of method would find the real root i.e. W_0(z)... This finds W_{-1}(z) local to the beam waist of an LG beam.
        
    end
    
    w = sqrt(-abs(paraxial_order )/2.*w); %Beam waist in normalized units
    
end
% knowing the angle, we approximate the waist scaling factor:
wscaling=1/tan(abs(anglez/180*pi));

%nmax is to estimate relevant parameters.
nmax = ceil( 2 * (pi * abs(w) + 4 + wscaling + paraxial_order));

ntheta=round(2*(itlim/200)^2*(lg_mode(1)+1)*(1/16^2*(nmax)^2+4*(nmax+2))+1);
delta=truncation_angle/ntheta/2;

[theta,phi]=meshgrid(linspace(delta,truncation_angle-delta,ntheta-1)/180*pi,0);

%LG MODES!!!
rw = (wscaling * w)^2 * tan(theta).^2 ;

L = laguerre(lg_mode(1),abs(lg_mode(2)),2*rw);
beam_envelope = sin(theta).*rw.^abs(lg_mode(2)/2) .* L .* exp(-rw + 1i*lg_mode(2)*phi + 1i*pi/2*(paraxial_order+1));

if debug
    figure(1010101)
    plot(theta*180/pi,abs(beam_envelope),'x');
end

Esmall=[beam_envelope.',beam_envelope.'];

%Now generate the modes:
m=[lg_mode(2)-1,lg_mode(2)+1];

if strcmp(class(polarisation),'char')
    m=lg_mode(2);
end

[nv,mv]=meshgrid([1:max([itlim,nmax])],m);

nn=nv(abs(mv)<=nv);
mm=mv(abs(mv)<=nv);

index=1;
a=zeros(length(nv(:)),1);
b=a;

%parameters for loop control
goooo=1;
ii=0;
currpwr=zeros(length(nv(:)),1);
while goooo
    ii=ii+1;
    mt=m(abs(m)<=ii);
    
    lmt=length(mt);
    
    [B,C,~]=vsh(ii,mt,theta,phi);
    for jj=1:lmt
        
        if abs(truncation_angle)<=90
            at = trapz(Esmall.*conj(C(:,jj+(lmt:lmt:end-lmt)))).';
            bt = trapz(Esmall.*conj(B(:,jj+(lmt:lmt:end-lmt)))).';
        else
            at = newtoncotes(Esmall.*conj(C(:,jj+(lmt:lmt:end-lmt)))).';
            bt = newtoncotes(Esmall.*conj(B(:,jj+(lmt:lmt:end-lmt)))).';
        end
        
        a(index)=[at(1)+(mt(jj)-lg_mode(2))*1i*at(2)]*1i^(ii+1)/sqrt(ii*(ii+1));
        b(index)=[bt(1)+(mt(jj)-lg_mode(2))*1i*bt(2)]*1i^(ii)/sqrt(ii*(ii+1));
        index=index+1;
        currpwr(ii)=sqrt(sum(abs(a).^2+abs(b).^2));
    end
    
    if ii>nmax+2%ii>=itlim%
        goooo=0;
        if and(mean((currpwr(ii-nmax+1:ii)-currpwr(ii-nmax:ii-1))./currpwr(ii-nmax+1:ii))>tol , ii<itlim)
            goooo=1;
        end
    end
    
end

if ii>=itlim
    warning('ott:bsc_lgmode_farfield:limit',['Iteration limit reached. Tolerance: ' num2str(tol) ', Current deviation: ' num2str(mean((currpwr(ii-nmax+1:ii)-currpwr(ii-nmax:ii-1))./currpwr(ii-nmax+1:ii))) '. Increase tolerance or iteration limit']);
end

if debug
    figure(1010102)
    plot(currpwr);
end

a=a(1:index-1);
b=b(1:index-1);
nn=nn(1:index-1);
mm=mm(1:index-1);

pplus=sqrt(sum(abs(a(mm==lg_mode(2)+1)).^2+abs(b(mm==lg_mode(2)+1)).^2));
pminus=sqrt(sum(abs(a(mm==lg_mode(2)-1)).^2+abs(b(mm==lg_mode(2)-1)).^2));

if pplus~=0
    a(mm==lg_mode(2)+1)=a(mm==lg_mode(2)+1)/pplus;
    b(mm==lg_mode(2)+1)=b(mm==lg_mode(2)+1)/pplus;
end
if pminus~=0
    a(mm==lg_mode(2)-1)=a(mm==lg_mode(2)-1)/pminus;
    b(mm==lg_mode(2)-1)=b(mm==lg_mode(2)-1)/pminus;
end

%Perform polarisation normalisation here.
if strcmp(class(polarisation),'char')
    switch polarisation
        case 'az'
            a=0*a;
        case 'ra'
            b=0*b;
    end
else
    
    T=-[1 1i;1 -1i]*((-1)^(lg_mode(2)))/sqrt(2); %Transformation matrix from jones vector to BSC
    c=T*polarisation(:)/sqrt(sum(abs(polarisation).^2));
    
    if debug
        c
    end
    
    a(mm==mt(1))=c(1)*a(mm==mt(1));
    a(mm==mt(2))=c(2)*a(mm==mt(2));
    
    b(mm==mt(1))=c(1)*b(mm==mt(1));
    b(mm==mt(2))=c(2)*b(mm==mt(2));
    
end

%output the modes of power.
pwr=sqrt(sum(abs(a).^2+abs(b).^2));
a=a/pwr;
b=b/pwr;

if nargout==2
    ci=combined_index(nn,mm);
    
    nn=sparse(max(nn)*(max(nn)+2),1);
    mm=nn;
    
    nn(ci)=a;
    mm(ci)=b;
end

return

