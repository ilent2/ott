function [structureoutput]=electromagnetic_field_xyz(kxyz,nm,ab,pq,cd,varargin)
%ELECTROMAGNETIC_FIELD_XYZ calculates the fields from beam vectors.
%
% S = ELECTROMAGNETIC_FIELD_XYZ(kr, nm, ab, pq, cd)
% calculates the field at points kr for beams ab, pq and cd with mode
% incides nm.
%
% kr is a vector of cartesian coordinates (times wavenumber).
% nm = [n;m] is a column vector of mode indices.
% ab = [a;b] is a column vector of full incident beam shape coefficients.
% pq = [p;q] is a column vector of full scattered beam shape coefficients.
% cd = [c;d] is a column vector of full internal beam shape coefficients.
%
% The beam vectors can be full or sparse column vectors.
%
% n,m are the mode indices, these can be in truncated form. The calculation
% will be quicker if a truncated n and m can be used.
%
% Any combination of beams can be calculated, beams can be omitted
% by replacing the corresponding beam with an empty array.
%
% The output is a structure which contains fields for each of the
% requested beams:
%     ab -> S.Eincident and S.Hincident
%     pq -> S.Escattered and S.Hscattered
%     cd -> S.Einternal and S.Hinternal
%
% S = ELECTROMAGNETIC_FIELD_XYZ(..., 'relativerefractiveindex', relidx)
% specifies the relative refractive index for internal field calculations.
%
% S = ELECTROMAGNETIC_FIELD_XYZ(..., 'tolerance', tol)
% may be used in future.
%
% S = ELECTROMAGNETIC_FIELD_XYZ(..., 'displacementfield') and
% S = ELECTROMAGNETIC_FIELD_XYZ(..., 'displacementfield', dis)
% may be used in future.
%
% NOTE: If internal fields are calculated only the theta and phi components
% of E are continuous at the boundary. Conversely, only the kr component of
% D is continuous at the boundary.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

verbose=0;

relindx=1;
tol=1e-8;
dfield=0;
if numel(varargin)>0
    for ii=1:length(varargin)
        switch class(varargin{ii})
            case 'cell'
                switch lower(varargin{ii}{1})
                    case 'relativerefractiveindex'
                        relindx=varargin{ii}{2};
                    case 'tolerance'
                        tol=varargin{ii}{2};
                    case 'displacementfield'
                        try
                            dfield=varargin{ii}{2};
                        catch
                            dfield=1;
                        end
                    otherwise
                      ott_warning('ott:electromagnetic_field_xyz:input', ...
                          ['Unrecognised input: ' varargin{ii} '.'])
                end
            otherwise
        end
    end
end

ott_warning('internal');

lengthnm=size(nm,1)/2;

n=nm(1:size(nm,1)/2,1);
m=nm(size(nm,1)/2+1:size(nm,1),1);

ci=combined_index(n,m);

try
    lengthab=size(ab,1)/2;
    ab=full(ab([ci;ci+lengthab]));
catch
    ab=[];
end

try
    lengthpq=size(pq,1)/2;
    pq=full(pq([ci;ci+lengthpq]));
catch
    pq=[];
end

try
    lengthcd=size(cd,1)/2;
    cd=full(cd([ci;ci+lengthcd]));
catch
    cd=[];
end

%calculate the space
lengthab=size(ab,1)/2;
lengthcd=size(cd,1)/2;
lengthpq=size(pq,1)/2;

[rv,tv,pv]=xyz2rtp(kxyz(:,1),kxyz(:,2),kxyz(:,3));

[r_new,~,indR]=unique(rv);
[theta_new,~,indTheta]=unique(tv);
[phi_new,~,indPhi]=unique(pv);

% r_new=r_new;
r_new(r_new==0)=1e-15;
%ab length can be bigger than cd or pq as it is the beam we start with.

%look for biggest and smallest elements in the matrix, kill elements
%less than tol of the max.

behaviour=0;
if lengthab~=0
    behaviour=behaviour+1;
    a=ab(1:lengthab,:);
    b=ab(lengthab+1:end,:);
end

if lengthpq~=0
    behaviour=behaviour+2;
    p=pq(1:lengthpq,:);
    q=pq(lengthpq+1:end,:);
end

if lengthcd~=0
    behaviour=behaviour+4;
    c=cd(1:lengthcd,:);
    d=cd(lengthcd+1:end,:);
end

if verbose
    behaviour
end

if behaviour==0
    ott_warning('external');
    error('no non-zero elements!')
end

E1 = zeros(size(kxyz));
E2 = E1;
E3 = E1;
H1 = E1;
H2 = E1;
H3 = E1;

un=unique(n);

%NumberOfLoops=length(calcvecels)
if verbose
    tic
end

switch behaviour
    case 1
        
        for nn = 1:max(un)
            if verbose
                nn
            end
            
            Nn = 1/sqrt(nn*(nn+1));
            vv=find(n==nn);
            
            if ~isempty(vv)
                
                [Y,Ytheta,Yphi] = spharm(nn,m(vv),theta_new,zeros(size(theta_new)));
                [jn,djn]=sbesselj(nn,r_new);
                
                [M,PHI]=meshgrid(m(vv),phi_new);
                
                expimphi=exp(1i*M.*PHI);%repmat(,[1,3]);
                
                jnU=jn(indR);
                djnU=djn(indR);
                kr=r_new(indR);
                
                for ii=1:length(vv)
                    index=vv(ii);%nn*(nn+1)+m(vv(ii))
                    
                    E1(:,1)=E1(:,1)+Nn*b(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E1(:,2)=E1(:,2)+Nn*(a(index)*jnU.*Yphi(indTheta,ii)+b(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E1(:,3)=E1(:,3)+Nn*(-a(index)*jnU.*Ytheta(indTheta,ii)+b(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    %if nargout>1
                    H1(:,1)=H1(:,1)+Nn*a(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H1(:,2)=H1(:,2)+Nn*(b(index)*jnU.*Yphi(indTheta,ii)+a(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H1(:,3)=H1(:,3)+Nn*(-b(index)*jnU.*Ytheta(indTheta,ii)+a(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    %end
                    
                end
            end
            
            
        end
        
    case 2
        
        for nn = 1:max(un)
            if verbose
                nn
            end
  
           Nn = 1/sqrt(nn*(nn+1));
           vv=find(n==nn);
            
            if ~isempty(vv)
                [Y,Ytheta,Yphi] = spharm(nn,m(vv),theta_new,zeros(size(theta_new)));
                
                [hn,dhn]=sbesselh1(nn,r_new);
                
                [M,PHI]=meshgrid(1i*m(vv),phi_new);
                
                expimphi=repmat(exp(M.*PHI),[1,3]);
                
                Nn = 1/sqrt(nn*(nn+1));
                
                
                hnU=hn(indR);
                dhnU=dhn(indR);
                
                kr=r_new(indR);
                
                for ii=1:length(vv)
                    index=vv(ii);%nn*(nn+1)+m(vv(ii))
                    
                    E2(:,1)=E2(:,1)+Nn*q(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E2(:,2)=E2(:,2)+Nn*(p(index)*hnU.*Yphi(indTheta,ii)+q(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E2(:,3)=E2(:,3)+Nn*(-p(index)*hnU.*Ytheta(indTheta,ii)+q(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    %                     if nargout>1
                    H2(:,1)=H2(:,1)+Nn*p(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H2(:,2)=H2(:,2)+Nn*(q(index)*hnU.*Yphi(indTheta,ii)+p(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H2(:,3)=H2(:,3)+Nn*(-q(index)*hnU.*Ytheta(indTheta,ii)+p(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    %                     end
                    
                end
            end
            
            
        end
        
    case 3
        
        for nn = 1:max(un)
            if verbose
                nn
            end
            Nn = 1/sqrt(nn*(nn+1));
            vv=find(n==nn);
            
            if ~isempty(vv)
                [Y,Ytheta,Yphi] = spharm(nn,m(vv),theta_new,zeros(size(theta_new)));
                [jn,djn]=sbesselj(nn,r_new);
                [hn,dhn]=sbesselh1(nn,r_new);
                
                [M,PHI]=meshgrid(1i*m(vv),phi_new);
                
                expimphi=repmat(exp(M.*PHI),[1,3]);
                                
                jnU=jn(indR);
                djnU=djn(indR);
                
                hnU=hn(indR);
                dhnU=dhn(indR);
                
                kr=r_new(indR);
                
                for ii=1:length(vv)
                    index=vv(ii);%nn*(nn+1)+m(vv(ii))
                    
                    E1(:,1)=E1(:,1)+Nn*b(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E1(:,2)=E1(:,2)+Nn*(a(index)*jnU.*Yphi(indTheta,ii)+b(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E1(:,3)=E1(:,3)+Nn*(-a(index)*jnU.*Ytheta(indTheta,ii)+b(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    E2(:,1)=E2(:,1)+Nn*q(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E2(:,2)=E2(:,2)+Nn*(p(index)*hnU.*Yphi(indTheta,ii)+q(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E2(:,3)=E2(:,3)+Nn*(-p(index)*hnU.*Ytheta(indTheta,ii)+q(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    %                     if nargout>1
                    H1(:,1)=H1(:,1)+Nn*a(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H1(:,2)=H1(:,2)+Nn*(b(index)*jnU.*Yphi(indTheta,ii)+a(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H1(:,3)=H1(:,3)+Nn*(-b(index)*jnU.*Ytheta(indTheta,ii)+a(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    H2(:,1)=H2(:,1)+Nn*p(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H2(:,2)=H2(:,2)+Nn*(q(index)*hnU.*Yphi(indTheta,ii)+p(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H2(:,3)=H2(:,3)+Nn*(-q(index)*hnU.*Ytheta(indTheta,ii)+p(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    %                     end
                    
                end
            end
            
            
        end
    case 4
        for nn = 1:max(un)
            if verbose
                nn
            end
            Nn = 1/sqrt(nn*(nn+1));
            vv=find(n==nn);
            
            if ~isempty(vv)
                [Y,Ytheta,Yphi] = spharm(nn,m(vv),theta_new,zeros(size(theta_new)));
                
                [jn,djn]=sbesselj(nn,r_new);
                [jnr,djnr]=sbesselj(nn,relindx*r_new); %relindx*
                [hn,dhn]=sbesselh1(nn,r_new);
                
                [M,PHI]=meshgrid(1i*m(vv),phi_new);
                
                expimphi=repmat(exp(M.*PHI),[1,3]);
                
                jnU=jn(indR);
                djnU=djn(indR);
                
                jnrU=jnr(indR);
                djnrU=djnr(indR);
                
                hnU=hn(indR);
                dhnU=dhn(indR);
                
                kr=r_new(indR);
                
                for ii=1:length(vv)
                    index=vv(ii);%nn*(nn+1)+m(vv(ii))
                    
                    E3(:,1)=E3(:,1)+Nn*d(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E3(:,2)=E3(:,2)+Nn*(c(index)*jnrU.*Yphi(indTheta,ii)+d(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E3(:,3)=E3(:,3)+Nn*(-c(index)*jnrU.*Ytheta(indTheta,ii)+d(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    H3(:,1)=H3(:,1)+Nn*c(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H3(:,2)=H3(:,2)+Nn*(d(index)*jnrU.*Yphi(indTheta,ii)+c(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H3(:,3)=H3(:,3)+Nn*(-d(index)*jnrU.*Ytheta(indTheta,ii)+c(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    %                     end
                    
                end
            end
            
            
        end
    case 5
        for nn = 1:max(un)
            if verbose
                nn
            end
            Nn = 1/sqrt(nn*(nn+1));
            vv=find(n==nn);
            
            if ~isempty(vv)
                [Y,Ytheta,Yphi] = spharm(nn,m(vv),theta_new,zeros(size(theta_new)));
                
                [jn,djn]=sbesselj(nn,r_new);
                [jnr,djnr]=sbesselj(nn,relindx*r_new); %relindx*
                
                [M,PHI]=meshgrid(1i*m(vv),phi_new);
                
                expimphi=repmat(exp(M.*PHI),[1,3]);
                
                jnU=jn(indR);
                djnU=djn(indR);
                
                jnrU=jnr(indR);
                djnrU=djnr(indR);
                
                kr=r_new(indR);
                
                for ii=1:length(vv)
                    index=vv(ii);%nn*(nn+1)+m(vv(ii))
                    
                    E1(:,1)=E1(:,1)+Nn*b(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E1(:,2)=E1(:,2)+Nn*(a(index)*jnU.*Yphi(indTheta,ii)+b(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E1(:,3)=E1(:,3)+Nn*(-a(index)*jnU.*Ytheta(indTheta,ii)+b(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    E3(:,1)=E3(:,1)+Nn*d(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E3(:,2)=E3(:,2)+Nn*(c(index)*jnrU.*Yphi(indTheta,ii)+d(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E3(:,3)=E3(:,3)+Nn*(-c(index)*jnrU.*Ytheta(indTheta,ii)+d(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    %                     if nargout>1
                    H1(:,1)=H1(:,1)+Nn*a(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H1(:,2)=H1(:,2)+Nn*(b(index)*jnU.*Yphi(indTheta,ii)+a(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H1(:,3)=H1(:,3)+Nn*(-b(index)*jnU.*Ytheta(indTheta,ii)+a(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    H3(:,1)=H3(:,1)+Nn*c(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H3(:,2)=H3(:,2)+Nn*(d(index)*jnrU.*Yphi(indTheta,ii)+c(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H3(:,3)=H3(:,3)+Nn*(-d(index)*jnrU.*Ytheta(indTheta,ii)+c(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    %                     end
                    
                end
            end
            
            
        end
    case 6
        for nn = 1:max(un)
            if verbose
                nn
            end
            Nn = 1/sqrt(nn*(nn+1));
            vv=find(n==nn);
            
            if ~isempty(vv)
                [Y,Ytheta,Yphi] = spharm(nn,m(vv),theta_new,zeros(size(theta_new)));
                
                [jnr,djnr]=sbesselj(nn,relindx*r_new); %relindx*
                [hn,dhn]=sbesselh1(nn,r_new);
                
                [M,PHI]=meshgrid(1i*m(vv),phi_new);
                
                expimphi=repmat(exp(M.*PHI),[1,3]);
                
                jnrU=jnr(indR);
                djnrU=djnr(indR);
                
                hnU=hn(indR);
                dhnU=dhn(indR);
                
                kr=r_new(indR);
                
                for ii=1:length(vv)
                    index=vv(ii);%nn*(nn+1)+m(vv(ii))
                    

                    E2(:,1)=E2(:,1)+Nn*q(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E2(:,2)=E2(:,2)+Nn*(p(index)*hnU.*Yphi(indTheta,ii)+q(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E2(:,3)=E2(:,3)+Nn*(-p(index)*hnU.*Ytheta(indTheta,ii)+q(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    E3(:,1)=E3(:,1)+Nn*d(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E3(:,2)=E3(:,2)+Nn*(c(index)*jnrU.*Yphi(indTheta,ii)+d(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E3(:,3)=E3(:,3)+Nn*(-c(index)*jnrU.*Ytheta(indTheta,ii)+d(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    %                     if nargout>1
            
                    H2(:,1)=H2(:,1)+Nn*p(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H2(:,2)=H2(:,2)+Nn*(q(index)*hnU.*Yphi(indTheta,ii)+p(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H2(:,3)=H2(:,3)+Nn*(-q(index)*hnU.*Ytheta(indTheta,ii)+p(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    H3(:,1)=H3(:,1)+Nn*c(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H3(:,2)=H3(:,2)+Nn*(d(index)*jnrU.*Yphi(indTheta,ii)+c(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H3(:,3)=H3(:,3)+Nn*(-d(index)*jnrU.*Ytheta(indTheta,ii)+c(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    %                     end
                    
                end
            end
            
            
        end
    case 7
        for nn = 1:max(un)
            if verbose
                nn
            end
            Nn = 1/sqrt(nn*(nn+1));
            vv=find(n==nn);
            
            if ~isempty(vv)
                [Y,Ytheta,Yphi] = spharm(nn,m(vv),theta_new,zeros(size(theta_new)));
                
                [jn,djn]=sbesselj(nn,r_new);
                [jnr,djnr]=sbesselj(nn,relindx*r_new); %relindx*
                [hn,dhn]=sbesselh1(nn,r_new);
                
                [M,PHI]=meshgrid(1i*m(vv),phi_new);
                
                expimphi=repmat(exp(M.*PHI),[1,3]);
                
                jnU=jn(indR);
                djnU=djn(indR);
                
                jnrU=jnr(indR);
                djnrU=djnr(indR);
                
                hnU=hn(indR);
                dhnU=dhn(indR);
                
                kr=r_new(indR);
                
                for ii=1:length(vv)
                    index=vv(ii);%nn*(nn+1)+m(vv(ii))
                    
                    E1(:,1)=E1(:,1)+Nn*b(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E1(:,2)=E1(:,2)+Nn*(a(index)*jnU.*Yphi(indTheta,ii)+b(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E1(:,3)=E1(:,3)+Nn*(-a(index)*jnU.*Ytheta(indTheta,ii)+b(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    E2(:,1)=E2(:,1)+Nn*q(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E2(:,2)=E2(:,2)+Nn*(p(index)*hnU.*Yphi(indTheta,ii)+q(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E2(:,3)=E2(:,3)+Nn*(-p(index)*hnU.*Ytheta(indTheta,ii)+q(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    E3(:,1)=E3(:,1)+Nn*d(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E3(:,2)=E3(:,2)+Nn*(c(index)*jnrU.*Yphi(indTheta,ii)+d(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E3(:,3)=E3(:,3)+Nn*(-c(index)*jnrU.*Ytheta(indTheta,ii)+d(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    %                     if nargout>1
                    H1(:,1)=H1(:,1)+Nn*a(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H1(:,2)=H1(:,2)+Nn*(b(index)*jnU.*Yphi(indTheta,ii)+a(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H1(:,3)=H1(:,3)+Nn*(-b(index)*jnU.*Ytheta(indTheta,ii)+a(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    H2(:,1)=H2(:,1)+Nn*p(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H2(:,2)=H2(:,2)+Nn*(q(index)*hnU.*Yphi(indTheta,ii)+p(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H2(:,3)=H2(:,3)+Nn*(-q(index)*hnU.*Ytheta(indTheta,ii)+p(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    H3(:,1)=H3(:,1)+Nn*c(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H3(:,2)=H3(:,2)+Nn*(d(index)*jnrU.*Yphi(indTheta,ii)+c(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H3(:,3)=H3(:,3)+Nn*(-d(index)*jnrU.*Ytheta(indTheta,ii)+c(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    %                     end
                    
                end
            end
            
            
        end
        
        
        
    otherwise
        ott_warning('external');
        error('ott:electromagnetic_field_xyz:notyetimplemented','This behaviour is not implemented.')
        
        
        
end
if verbose
    toc
end

 H1=-1i*H1; %LOOK HERE TO FIX STUFF
 H2=-1i*H2; %LOOK HERE TO FIX STUFF
 H3=-1i*H3; %LOOK HERE TO FIX STUFF


%res flipped because it's a meshgrid
if behaviour==7||behaviour==4||behaviour==5||behaviour==6
    
    [Exv3,Eyv3,Ezv3] = rtpv2xyzv(E3(:,1),E3(:,2),E3(:,3),rv,tv,pv);
    Ex3=Exv3;
    index= isnan(Ex3);
    Ex3(index)=0;
    Ey3=Eyv3;
    index= isnan(Ey3);
    Ey3(index)=0;
    Ez3=Ezv3;
    index= isnan(Ez3);
    Ez3(index)=0;
    
    [Hxv3,Hyv3,Hzv3] = rtpv2xyzv(H3(:,1),H3(:,2),H3(:,3),rv,tv,pv);
    Hx3=Hxv3;
    index= isnan(Hx3);
    Hx3(index)=0;
    Hy3=Hyv3;
    index= isnan(Hy3);
    Hx3(index)=0;
    Hz3=Hzv3;
    index= isnan(Hz3);
    Hx3(index)=0;
    
    structureoutput.Einternal=[squeeze(Ex3),squeeze(Ey3),squeeze(Ez3)];
    structureoutput.Hinternal=[squeeze(Hx3),squeeze(Hy3),squeeze(Hz3)];
    
    ott_warning('ott:electromagnetic_field_xyz:internalfield', ...
        'Must scale grid for internal fields by relative refractive index.')
end

if behaviour==2|behaviour==3|behaviour==6|behaviour==7
    
    [Exv2,Eyv2,Ezv2] = rtpv2xyzv(E2(:,1),E2(:,2),E2(:,3),rv,tv,pv);
    Ex2=Exv2;
    index= isnan(Ex2);
    Ex2(index)=0;
    Ey2=Eyv2;
    index= isnan(Ey2);
    Ey2(index)=0;
    Ez2=Ezv2;
    index= isnan(Ez2);
    Ez2(index)=0;
    
    [Hxv2,Hyv2,Hzv2] = rtpv2xyzv(H2(:,1),H2(:,2),H2(:,3),rv,tv,pv);
    Hx2=Hxv2;
    index= isnan(Hx2);
    Hx2(index)=0;
    Hy2=Hyv2;
    index= isnan(Hy2);
    Hx2(index)=0;
    Hz2=Hzv2;
    index= isnan(Hz2);
    Hx2(index)=0;
    
    structureoutput.Escattered=[squeeze(Ex2),squeeze(Ey2),squeeze(Ez2)];
    structureoutput.Hscattered=[squeeze(Hx2),squeeze(Hy2),squeeze(Hz2)];
    
end
if behaviour==1|behaviour==3|behaviour==5|behaviour==7
    
    [Exv1,Eyv1,Ezv1] = rtpv2xyzv(E1(:,1),E1(:,2),E1(:,3),rv,tv,pv);
    Ex1=Exv1;
    index= isnan(Ex1);
    Ex1(index)=0;
    Ey1=Eyv1;
    index= isnan(Ey1);
    Ey1(index)=0;
    Ez1=Ezv1;
    index= isnan(Ez1);
    Ez1(index)=0;
    
    [Hxv1,Hyv1,Hzv1] = rtpv2xyzv(H1(:,1)*1,1*H1(:,2),1*H1(:,3),rv,tv,pv);
    Hx1=Hxv1;
    index= isnan(Hx1);
    Hx1(index)=0;
    Hy1=Hyv1;
    index= isnan(Hy1);
    Hx1(index)=0;
    Hz1=Hzv1;
    index= isnan(Hz1);
    Hx1(index)=0;
    
    structureoutput.Eincident=[squeeze(Ex1),squeeze(Ey1),squeeze(Ez1)];
    structureoutput.Hincident=[squeeze(Hx1),squeeze(Hy1),squeeze(Hz1)];
    
end

ott_warning('external');

return
