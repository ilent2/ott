function [nn,mm,a,b] = bsc_hgmode_farfield( anglez,hg_mode,varargin )
% bsc_hgmode_farfield.m
% Uses a line integral method of LG beams to find the
% spherical harmonic expansion of a Hermite-Gaussian laser beam.
%
% Usage:
% [n,m,a,b] = bsc_hgmode_farfield( angle,hg_mode );
%
% calculates the expansion coefficients for a 'hg_mode' ([m,n]) x-polarised
% beam converging with a cone of 'angle' (in degrees)
%
% OR
%
% [n,m,a,b] = bsc_hgmode_farfield( angle,hg_mode,(optional))
%
% (optional) inputs are not limited to a single entry or order and include:
%
% truncation_angle --- angle of hard edge aperture (degrees) default: 90
% polarisation --- either a Jones vetor {[E_x,Ey] | 'azimuthal' |
% 'radial'}. default: [1,0]
% NYI!!!: mapping --- {'direct' | {'projection'}}. 'direct' uses non-paraxial
% mapping. 'projection' uses the symmetry maintinging
% paraxial-to-non-paraxial convention.
%
% PACKAGE INFO

w0_scaling='free';
nmax=200;

%setting up
if nargin>2
    for ii=1:length(varargin)
        switch class(varargin{ii})
            case {'char'}
                switch lower(varargin{ii})
                    case 'fixed'
                        w0_scaling='fixed';
                end
        end
    end
end

paraxial_order=sum(hg_mode);

[modeweights,lgindex,hgindex]=genLG2HG(paraxial_order);

[m_,n_]=hglookup(paraxial_order,hgindex);

[row,col]=find(m_==hg_mode(1),1);

% pick angle based on max l mode...
if strcmp(lower(w0_scaling),'free')
    w = 1.; %Beam waist in normalized units.
    
    if (order ~= 0)
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
    
    anglez=atan2(tan(anglez/180*pi),w)*180/pi;
end

att=sparse(nmax*(nmax+2),1);
btt=att;

for ii=1:paraxial_order+1
    
    [p,l]=lglookup(paraxial_order,lgindex(row,ii));
    
    [n,m,at,bt]=bsc_lgmode_farfield(anglez,[p,l],'fixed',varargin{:});
    
    ci=combined_index(n,m);
    
    att(ci)=att(ci)+modeweights(row,ii)*at;
    btt(ci)=btt(ci)+modeweights(row,ii)*bt;
    
end

if nargout==2
    
    inda=find(att,1,'last');
    indb=find(btt,1,'last');
    
    nmax=floor(sqrt(max([inda,indb])));
    
    nn=att(1:nmax*(nmax+2));
    mm=btt(1:nmax*(nmax+2));
    
else
    ci=unique([find(att);find(btt)]);
    
    [nn,mm]=combined_index(ci(:));
    
    a=full(att(ci));
    b=full(btt(ci));
end