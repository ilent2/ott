function [structureoutput, data]=electromagnetic_field_xyz(kxyz,nm, ...
    ibeam,obeam,rbeam,varargin)
%ELECTROMAGNETIC_FIELD_XYZ calculates the fields from beam vectors.
%
% S = ELECTROMAGNETIC_FIELD_XYZ(kr, nm, ibeam, obeam, rbeam)
% calculates the field at points kr for beams with mode incides nm.
%
% [S, data] = ELECTROMAGNETIC_FIELD_XYZ(..., 'saveData', true) saves
% data that can be used for repeated calculations.
%
% kr is a vector of cartesian coordinates (times wavenumber).
% ibeam is the incident beam, obeam is the scattered beam, rbeam
% is the internal beam.
%
% n,m are the mode indices, these can be in truncated form. The calculation
% will be quicker if a truncated n and m can be used.
%
% Any combination of beams can be calculated, beams can be omitted
% by replacing the corresponding beam with an empty array.
%
% The output is a structure which contains fields for each of the
% requested beams:
%     ibeam -> S.Eincident and S.Hincident
%     obeam -> S.Escattered and S.Hscattered
%     rbeam -> S.Einternal and S.Hinternal
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
% Optional named arguments:
%     'calcE'   bool   calculate E field (default: true)
%     'calcH'   bool   calculate H field (default: true)
%
% NOTE: If internal fields are calculated only the theta and phi components
% of E are continuous at the boundary. Conversely, only the kr component of
% D is continuous at the boundary.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

import ott.utils.*

ip = inputParser;
ip.addParameter('relativerefractiveindex', 1);
ip.addParameter('tolerance', 1e-8);
ip.addParameter('displacementfield', 0);
ip.addParameter('verbose', false);
ip.addParameter('calcE', true);
ip.addParameter('calcH', true);
ip.addParameter('saveData', nargout == 2);
ip.addParameter('data', []);
ip.parse(varargin{:});

relindx = ip.Results.relativerefractiveindex;
% tol = ip.Results.tolerance;
% dfield = ip.Results.displacementfield;
verbose = ip.Results.verbose;

ott.warning('internal');

n=nm(1:size(nm,1)/2,1);
m=nm(size(nm,1)/2+1:size(nm,1),1);

ci=combined_index(n,m);

behaviour=0;

if ~isempty(ibeam)
  ibeam.Nmax = max(n, ibeam.Nmax);
  [a,b] = ibeam.getCoefficients(ci);
  behaviour=behaviour+1;
end

if ~isempty(obeam)
  obeam.Nmax = max(n, obeam.Nmax);
  [p,q] = obeam.getCoefficients(ci);
  behaviour=behaviour+2;
end

if ~isempty(rbeam)
  rbeam.Nmax = max(n, rbeam.Nmax);
  [c,d] = rbeam.getCoefficients(ci);
  behaviour=behaviour+4;
end

if behaviour==0
    ott.warning('external');
    error('no non-zero elements!')
end

[rv,tv,pv]=xyz2rtp(kxyz(:,1),kxyz(:,2),kxyz(:,3));

[r_new,~,indR]=unique(rv);
[theta_new,~,indTheta]=unique(tv);
[phi_new,~,indPhi]=unique(pv);

% r_new=r_new;
r_new(r_new==0)=1e-15;
%ab length can be bigger than cd or pq as it is the beam we start with.

%look for biggest and smallest elements in the matrix, kill elements
%less than tol of the max.

if verbose
    disp(['emfieldxyz behaviour:', num2str(behaviour)]);
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
      
% Alocate memory for output data
data = [];
if ip.Results.saveData
  data = zeros(numel(indTheta), 0);
end

% Start a counter for accessing the data
if ~isempty(ip.Results.data)
  dataCount = 0;
end

switch behaviour
    case 1
        
        for nn = 1:max(un)
            if verbose
                disp(['emfieldxyz nn:', num2str(nn)]);
            end
            
            Nn = 1/sqrt(nn*(nn+1));
            vv=find(n==nn);
            
            if ~isempty(vv)
                
                [Y,Ytheta,Yphi] = spharm(nn,m(vv),theta_new,zeros(size(theta_new)));
                [jn,djn]=sbesselj(nn,r_new);
                
                [M,PHI]=meshgrid(m(vv),phi_new);
                
                expimphi=exp(1i*M.*PHI);%repmat(,[1,3]);
                
                % Create full Y, Ytheta, Yphi, expimphi matrices (opt, R2018a)
                Yf = Y(indTheta, :);
                Ythetaf = Ytheta(indTheta, :);
                Yphif = Yphi(indTheta, :);
                expimphif = expimphi(indPhi, :);
                
                jnU=jn(indR);
                djnU=djn(indR);
                kr=r_new(indR);
                
                % Now we use full matrices, we can use matmul (opt, R2018a)
                % TODO: Repeat this optimisation elsewhere
                if ip.Results.calcE
                  aidx = full(a(vv));
                  bidx = full(b(vv));
                  E1(:,1)=E1(:,1)+Nn*nn*(nn+1)./kr.*jnU.*(Yf.*expimphif)*bidx(:);
                  E1(:,2)=E1(:,2)+Nn*((jnU(:).*Yphif.*expimphif)*aidx(:)+(djnU(:).*Ythetaf.*expimphif)*bidx(:));
                  E1(:,3)=E1(:,3)+Nn*((-jnU(:).*Ythetaf.*expimphif)*aidx(:)+(djnU(:).*Yphif.*expimphif)*bidx(:));
                end
                
                for ii=1:length(vv)
                    index=vv(ii);%nn*(nn+1)+m(vv(ii))
                    
                    if ip.Results.calcE
%                     E1(:,1)=E1(:,1)+Nn*b(index)*nn*(nn+1)./kr.*jnU.*Yf(:,ii).*expimphif(:,ii);
%                     E1(:,2)=E1(:,2)+Nn*(a(index)*jnU.*Yphif(:,ii)+b(index)*djnU.*Ythetaf(:,ii)).*expimphif(:,ii);
%                     E1(:,3)=E1(:,3)+Nn*(-a(index)*jnU.*Ythetaf(:,ii)+b(index)*djnU.*Yphif(:,ii)).*expimphif(:,ii);
                    end
                    
                    if ip.Results.calcH
                    H1(:,1)=H1(:,1)+Nn*a(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H1(:,2)=H1(:,2)+Nn*(b(index)*jnU.*Yphi(indTheta,ii)+a(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H1(:,3)=H1(:,3)+Nn*(-b(index)*jnU.*Ytheta(indTheta,ii)+a(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    end
                    
                end
            end
            
            
        end
        
    case 2
        
        for nn = 1:max(un)
            if verbose
                disp(['emfieldxyz nn:', num2str(nn)]);
            end
  
%            Nn = 1/sqrt(nn*(nn+1));
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
                    
                    if ip.Results.calcE
                    E2(:,1)=E2(:,1)+Nn*q(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E2(:,2)=E2(:,2)+Nn*(p(index)*hnU.*Yphi(indTheta,ii)+q(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E2(:,3)=E2(:,3)+Nn*(-p(index)*hnU.*Ytheta(indTheta,ii)+q(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    end
                    
                    if ip.Results.calcH
                    H2(:,1)=H2(:,1)+Nn*p(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H2(:,2)=H2(:,2)+Nn*(q(index)*hnU.*Yphi(indTheta,ii)+p(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H2(:,3)=H2(:,3)+Nn*(-q(index)*hnU.*Ytheta(indTheta,ii)+p(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    end
                    
                end
            end
            
            
        end
        
    case 3
        
        for nn = 1:max(un)
            if verbose
                disp(['emfieldxyz nn:', num2str(nn)]);
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
                    
                    if ip.Results.calcE
                    E1(:,1)=E1(:,1)+Nn*b(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E1(:,2)=E1(:,2)+Nn*(a(index)*jnU.*Yphi(indTheta,ii)+b(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E1(:,3)=E1(:,3)+Nn*(-a(index)*jnU.*Ytheta(indTheta,ii)+b(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    E2(:,1)=E2(:,1)+Nn*q(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E2(:,2)=E2(:,2)+Nn*(p(index)*hnU.*Yphi(indTheta,ii)+q(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E2(:,3)=E2(:,3)+Nn*(-p(index)*hnU.*Ytheta(indTheta,ii)+q(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    end
                    
                    if ip.Results.calcH
                    H1(:,1)=H1(:,1)+Nn*a(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H1(:,2)=H1(:,2)+Nn*(b(index)*jnU.*Yphi(indTheta,ii)+a(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H1(:,3)=H1(:,3)+Nn*(-b(index)*jnU.*Ytheta(indTheta,ii)+a(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    H2(:,1)=H2(:,1)+Nn*p(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H2(:,2)=H2(:,2)+Nn*(q(index)*hnU.*Yphi(indTheta,ii)+p(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H2(:,3)=H2(:,3)+Nn*(-q(index)*hnU.*Ytheta(indTheta,ii)+p(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    end
                    
                end
            end
            
            
        end
    case 4
      
        for nn = 1:max(un)
            if verbose
                disp(['emfieldxyz nn:', num2str(nn)]);
            end
            Nn = 1/sqrt(nn*(nn+1));
            vv=find(n==nn);
            
            if ~isempty(vv)
              
              kr=r_new(indR);
              
              if isempty(ip.Results.data)
                
                [Y,Ytheta,Yphi] = spharm(nn,m(vv),theta_new,zeros(size(theta_new)));
                
                [jnr,djnr]=sbesselj(nn,relindx*r_new);
                
                [M,PHI]=meshgrid(1i*m(vv),phi_new);
                
                expimphi=exp(M.*PHI);
                
                jnrU=jnr(indR);
                djnrU=djnr(indR);
                
                % Create full Y, Ytheta, Yphi, expimphi matrices (opt, R2018a)
                expimphif = expimphi(indPhi, :);
                YExpf = Y(indTheta, :).*expimphif;
                YthetaExpf = Ytheta(indTheta, :).*expimphif;
                YphiExpf = Yphi(indTheta, :).*expimphif;
                
                % Save the data if requested
                if ip.Results.saveData
                  data(:, end+1) = jnrU;
                  data(:, end+1) = djnrU;
                  data(:, end+(1:size(Ytheta, 2))) = YExpf;
                  data(:, end+(1:size(Ytheta, 2))) = YthetaExpf;
                  data(:, end+(1:size(Ytheta, 2))) = YphiExpf;
                end
                
              else
                
                % Load the data if present
                jnrU = ip.Results.data(:, dataCount+1);
                dataCount = dataCount + 1;
                djnrU = ip.Results.data(:, dataCount+1);
                dataCount = dataCount + 1;
                YExpf = ip.Results.data(:, dataCount+(1:length(vv)));
                dataCount = dataCount + length(vv);
                YthetaExpf = ip.Results.data(:, dataCount+(1:length(vv)));
                dataCount = dataCount + length(vv);
                YphiExpf = ip.Results.data(:, dataCount+(1:length(vv)));
                dataCount = dataCount + length(vv);
                
              end
              
              cidx = full(c(vv));
              didx = full(d(vv));

              % Now we use full matrices, we can use matmul (opt, R2018a)
              if ip.Results.calcE
                E3(:,1)=E3(:,1)+Nn*nn*(nn+1)./kr./relindx.*jnrU.*(YExpf)*didx(:);
                E3(:,2)=E3(:,2)+Nn*((jnrU(:).*YphiExpf)*cidx(:)+(djnrU(:).*YthetaExpf)*didx(:));
                E3(:,3)=E3(:,3)+Nn*((-jnrU(:).*YthetaExpf)*cidx(:)+(djnrU(:).*YphiExpf)*didx(:));
              end
                
              if ip.Results.calcH
                H3(:,1)=H3(:,1)+Nn*nn*(nn+1)./kr./relindx.*jnrU.*YExpf*cidx(:);
                H3(:,2)=H3(:,2)+Nn*((jnrU(:).*YphiExpf)*didx(:)+(djnrU(:).*YthetaExpf)*cidx(:));
                H3(:,3)=H3(:,3)+Nn*((-jnrU(:).*YthetaExpf)*didx(:)+(djnrU(:).*YphiExpf)*cidx(:));
              end
              
            end
            
            
        end
    case 5
        for nn = 1:max(un)
            if verbose
                disp(['emfieldxyz nn:', num2str(nn)]);
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
                    
                    if ip.Results.calcE
                    E1(:,1)=E1(:,1)+Nn*b(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E1(:,2)=E1(:,2)+Nn*(a(index)*jnU.*Yphi(indTheta,ii)+b(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E1(:,3)=E1(:,3)+Nn*(-a(index)*jnU.*Ytheta(indTheta,ii)+b(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    E3(:,1)=E3(:,1)+Nn*d(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E3(:,2)=E3(:,2)+Nn*(c(index)*jnrU.*Yphi(indTheta,ii)+d(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E3(:,3)=E3(:,3)+Nn*(-c(index)*jnrU.*Ytheta(indTheta,ii)+d(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    end
                    
                    if ip.Results.calcH
                    H1(:,1)=H1(:,1)+Nn*a(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H1(:,2)=H1(:,2)+Nn*(b(index)*jnU.*Yphi(indTheta,ii)+a(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H1(:,3)=H1(:,3)+Nn*(-b(index)*jnU.*Ytheta(indTheta,ii)+a(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    H3(:,1)=H3(:,1)+Nn*c(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H3(:,2)=H3(:,2)+Nn*(d(index)*jnrU.*Yphi(indTheta,ii)+c(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H3(:,3)=H3(:,3)+Nn*(-d(index)*jnrU.*Ytheta(indTheta,ii)+c(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    end
                    
                end
            end
            
            
        end
    case 6
        for nn = 1:max(un)
            if verbose
                disp(['emfieldxyz nn:', num2str(nn)]);
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
                    

                    if ip.Results.calcE
                    E2(:,1)=E2(:,1)+Nn*q(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E2(:,2)=E2(:,2)+Nn*(p(index)*hnU.*Yphi(indTheta,ii)+q(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E2(:,3)=E2(:,3)+Nn*(-p(index)*hnU.*Ytheta(indTheta,ii)+q(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    E3(:,1)=E3(:,1)+Nn*d(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E3(:,2)=E3(:,2)+Nn*(c(index)*jnrU.*Yphi(indTheta,ii)+d(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E3(:,3)=E3(:,3)+Nn*(-c(index)*jnrU.*Ytheta(indTheta,ii)+d(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    end
                    
                    if ip.Results.calcH
                    H2(:,1)=H2(:,1)+Nn*p(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H2(:,2)=H2(:,2)+Nn*(q(index)*hnU.*Yphi(indTheta,ii)+p(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H2(:,3)=H2(:,3)+Nn*(-q(index)*hnU.*Ytheta(indTheta,ii)+p(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    H3(:,1)=H3(:,1)+Nn*c(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H3(:,2)=H3(:,2)+Nn*(d(index)*jnrU.*Yphi(indTheta,ii)+c(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H3(:,3)=H3(:,3)+Nn*(-d(index)*jnrU.*Ytheta(indTheta,ii)+c(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    end
                    
                end
            end
            
            
        end
    case 7
        for nn = 1:max(un)
            if verbose
                disp(['emfieldxyz nn:', num2str(nn)]);
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
                    
                    if ip.Results.calcE
                    E1(:,1)=E1(:,1)+Nn*b(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E1(:,2)=E1(:,2)+Nn*(a(index)*jnU.*Yphi(indTheta,ii)+b(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E1(:,3)=E1(:,3)+Nn*(-a(index)*jnU.*Ytheta(indTheta,ii)+b(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    E2(:,1)=E2(:,1)+Nn*q(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E2(:,2)=E2(:,2)+Nn*(p(index)*hnU.*Yphi(indTheta,ii)+q(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E2(:,3)=E2(:,3)+Nn*(-p(index)*hnU.*Ytheta(indTheta,ii)+q(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    E3(:,1)=E3(:,1)+Nn*d(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    E3(:,2)=E3(:,2)+Nn*(c(index)*jnrU.*Yphi(indTheta,ii)+d(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    E3(:,3)=E3(:,3)+Nn*(-c(index)*jnrU.*Ytheta(indTheta,ii)+d(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    end
                    
                    if ip.Results.calcH
                    H1(:,1)=H1(:,1)+Nn*a(index)*nn*(nn+1)./kr.*jnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H1(:,2)=H1(:,2)+Nn*(b(index)*jnU.*Yphi(indTheta,ii)+a(index)*djnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H1(:,3)=H1(:,3)+Nn*(-b(index)*jnU.*Ytheta(indTheta,ii)+a(index)*djnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    H2(:,1)=H2(:,1)+Nn*p(index)*nn*(nn+1)./kr.*hnU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H2(:,2)=H2(:,2)+Nn*(q(index)*hnU.*Yphi(indTheta,ii)+p(index)*dhnU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H2(:,3)=H2(:,3)+Nn*(-q(index)*hnU.*Ytheta(indTheta,ii)+p(index)*dhnU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    
                    H3(:,1)=H3(:,1)+Nn*c(index)*nn*(nn+1)./kr/relindx.*jnrU.*Y(indTheta,ii).*expimphi(indPhi,ii);
                    H3(:,2)=H3(:,2)+Nn*(d(index)*jnrU.*Yphi(indTheta,ii)+c(index)*djnrU.*Ytheta(indTheta,ii)).*expimphi(indPhi,ii);
                    H3(:,3)=H3(:,3)+Nn*(-d(index)*jnrU.*Ytheta(indTheta,ii)+c(index)*djnrU.*Yphi(indTheta,ii)).*expimphi(indPhi,ii);
                    end
                    
                end
            end
            
            
        end
        
        
        
    otherwise
        ott.warning('external');
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
    
    ott.warning('ott:electromagnetic_field_xyz:internalfield', ...
        'Must scale grid for internal fields by relative refractive index.')
end

if behaviour==2||behaviour==3||behaviour==6||behaviour==7
    
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
if behaviour==1||behaviour==3||behaviour==5||behaviour==7
    
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

ott.warning('external');
