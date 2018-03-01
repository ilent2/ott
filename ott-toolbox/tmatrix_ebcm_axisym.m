function T = tmatrix_ebcm_axisym(Nmax,k_medium,k_particle,rho,z);
% tmatrix_ebcm_axisym.m : extendended boundary condition method for
%                         axisymmetric particles.
%
% Usage:
% T = tmatrix_ebcm_axisym(Nmax,k_medium,k_particle,rho,z)
% rho is the cylindrical r coodinate of a corner or feature.
% z is the cylindrical z coordinate of a corner or feature.
%
% PACKAGE INFO

verbose=0;
% %%%% TEST VALUES
% verbose=1;
% Nmax=10;
% k_medium=2*pi*1;
% k_particle=2*pi*1.18;
% rho=[0,.25,.25,0];
% z=[.25,.25,-.25,-.25];
% % Tmie=tmatrix_mie(Nmax,k_medium,k_particle,.5);
% % rho=.5*sin(linspace(0,pi,13313));
% % z=.5*cos(linspace(0,pi,13313));
% %%%%

%generate the surface points and normals in spherical coordinates.
[rtp,n,ds]=axisym_boundarypoints(Nmax,rho,z);

%%%% set up containters
Y=cell(Nmax,1);
Ytheta=Y;
Yphi=Y;

Nn=zeros(Nmax,1);

kr=k_particle*rtp(:,1);
ikr=1./kr;
kr_=k_medium*rtp(:,1);
ikr_=1./kr_;

jkr=cell(Nmax,1);
jkr_=jkr;
hkr_=jkr_;

djkr=jkr;
djkr_=djkr;
dhkr_=djkr_;

Nn=1./sqrt([1:Nmax].*([1:Nmax]+1));

Nm=repmat(Nn,[Nmax,1]);
Nm=Nm.'.*Nm;

%make all the spharms and bessel functions NOW!!!!
for ii=1:Nmax
    
    [Y{ii},Ytheta{ii},Yphi{ii}] = spharm(ii,rtp(:,2),rtp(:,3));
    
    jkr{ii}=sbesselj(ii,kr);
    jkr_{ii}=sbesselj(ii,kr_);
    hkr_{ii}=sbesselh1(ii,kr_);
    
    djkr{ii}=sbesselj(ii-1,kr)-ii*jkr{ii}./kr;
    djkr_{ii}=sbesselj(ii-1,kr_)-ii*jkr_{ii}./kr_;
    dhkr_{ii}=sbesselh1(ii-1,kr_)-ii*hkr_{ii}./kr_;
    
end

lengthJs=Nmax^2*(Nmax+2)-sum([1:Nmax-1].*([1:Nmax-1]+1));
J11=zeros(lengthJs,1);
J12=J11;
J21=J11;
J22=J11;

RgJ11=J11;
RgJ12=J11;
RgJ21=J11;
RgJ22=J11;

i_indexes=J11;
j_indexes=J11;

Qv=zeros(lengthJs*4,1);
RgQv=zeros(lengthJs*4,1);

jcount=0;

if verbose
    tic;
end

fillnum=(2*Nmax+1);
rM=repmat(n(:,1),[1,fillnum]);
thetaM=repmat(n(:,2),[1,fillnum]);
iKr=repmat(ikr,[1,fillnum]);
iKr_=repmat(ikr_,[1,fillnum]);
dS=repmat(ds,[1,fillnum]);

pi2=2*pi;
%%%% end setup containers

for jj=1:Nmax
    %remember that when phi = 0, conj(Yphi)=-Yphi, conj(Ytheta)=Ytheta.
    
    for kk=1:Nmax
        p=min(kk,jj);
        
        indexzk=kk+1+[-p:p];
        indexzj=jj+1+[-p:p];
        
        els=[1:length(indexzk)];
        
        jindx=jcount+els;
        
        i_indexes(jindx)=kk*(kk+1)+[-p:p];
        j_indexes(jindx)=jj*(jj+1)+[-p:p];
        
	%Calculate the cross products of the spherical harmonics
        YtYp=-Ytheta{jj}(:,indexzj).*Yphi{kk}(:,indexzk);
        YpYt=Yphi{jj}(:,indexzj).*Ytheta{kk}(:,indexzk);
        
        YYt=Y{jj}(:,indexzj).*Ytheta{kk}(:,indexzk);
        YtY=Ytheta{jj}(:,indexzj).*Y{kk}(:,indexzk);
        
        YYp=-Y{jj}(:,indexzj).*Yphi{kk}(:,indexzk);
        YpY=Yphi{jj}(:,indexzj).*Y{kk}(:,indexzk);
        
        YpYp=-Yphi{jj}(:,indexzj).*Yphi{kk}(:,indexzk);
        YtYt=Ytheta{jj}(:,indexzj).*Ytheta{kk}(:,indexzk);
        
	%calculate the cross products of spherical bessel functions.
        jh_=repmat(jkr{jj}.*conj(hkr_{kk}),[1,fillnum]);
        jj_=repmat(jkr{jj}.*jkr_{kk},[1,fillnum]);
        
        jdh_=repmat(jkr{jj}.*conj(dhkr_{kk}),[1,fillnum]);
        djh_=repmat(djkr{jj}.*conj(hkr_{kk}),[1,fillnum]);
        
        jdj_=repmat(jkr{jj}.*djkr_{kk},[1,fillnum]);
        djj_=repmat(djkr{jj}.*jkr_{kk},[1,fillnum]);
        
        djdh_=repmat(djkr{jj}.*conj(dhkr_{kk}),[1,fillnum]);
        djdj_=repmat(djkr{jj}.*djkr_{kk},[1,fillnum]);

	%perform the cross product followed by dotting the normal vector and summing over the area elements.
        J11(jindx)=pi2*Nm(jj,kk)*sum(dS(:,els).*rM(:,els).*jh_(:,els).*(YtYp-YpYt),1);
        J12(jindx)=pi2*Nm(jj,kk)*sum(dS(:,els).*(rM(:,els).*jdh_(:,els).*(YpYp+YtYt) ...
				-kk*(kk+1)*iKr_(:,els).*thetaM(:,els).*jh_(:,els).*YtY),1);
        J21(jindx)=pi2*Nm(jj,kk)*sum(dS(:,els).*(-rM(:,els).*djh_(:,els).*(YtYt+YpYp) ... 
				+jj*(jj+1)*iKr(:,els).*thetaM(:,els).*jh_(:,els).*YYt),1);
        J22(jindx)=pi2*Nm(jj,kk)*sum(dS(:,els).*(rM(:,els).*djdh_(:,els).*(YtYp-YpYt) ... 
				-thetaM(:,els).*(jj*(jj+1)*iKr(:,els).*jdh_(:,els).*YYp-kk*(kk+1)*iKr_(:,els).*djh_(:,els).*YpY)),1);
        
        RgJ11(jindx)=pi2*Nm(jj,kk)*sum(dS(:,els).*rM(:,els).*jj_(:,els).*(YtYp-YpYt),1);
        RgJ12(jindx)=pi2*Nm(jj,kk)*sum(dS(:,els).*(rM(:,els).*jdj_(:,els).*(YpYp+YtYt) ... 
				-kk*(kk+1)*iKr_(:,els).*thetaM(:,els).*jj_(:,els).*YtY),1);
        RgJ21(jindx)=pi2*Nm(jj,kk)*sum(dS(:,els).*(-rM(:,els).*djj_(:,els).*(YtYt+YpYp) ... 
				+jj*(jj+1)*iKr(:,els).*thetaM(:,els).*jj_(:,els).*YYt),1);
        RgJ22(jindx)=pi2*Nm(jj,kk)*sum(dS(:,els).*(rM(:,els).*djdj_(:,els).*(YtYp-YpYt) ... 
				-thetaM(:,els).*(jj*(jj+1)*iKr(:,els).*jdj_(:,els).*YYp-kk*(kk+1)*iKr_(:,els).*djj_(:,els).*YpY)),1);
        
        jcount=jcount+length(indexzj);
        
    end
    if verbose
        disp(['n=' num2str(jj) ' at: ' num2str(toc) ' seconds.'])
    end
    
end

%Set up the Q and RgQ for sparse input:
Qv(1:lengthJs)=-1i*k_medium*(k_particle*J21+k_medium*J12);
Qv(lengthJs+1:lengthJs*2)=-1i*k_medium*(k_particle*J22+k_medium*J11);
Qv(2*lengthJs+1:lengthJs*3)=-1i*k_medium*(k_particle*J11+k_medium*J22); 
Qv(3*lengthJs+1:lengthJs*4)=-1i*k_medium*(k_particle*J12+k_medium*J21);

RgQv(1:lengthJs)=-1i*k_medium*(k_particle*RgJ21+k_medium*RgJ12);
RgQv(lengthJs+1:lengthJs*2)=-1i*k_medium*(k_particle*RgJ22+k_medium*RgJ11);
RgQv(2*lengthJs+1:lengthJs*3)=-1i*k_medium*(k_particle*RgJ11+k_medium*RgJ22);
RgQv(3*lengthJs+1:lengthJs*4)=-1i*k_medium*(k_particle*RgJ12+k_medium*RgJ21);

%Convert the vectors of Q and RgQ into sparse vector form:
Q=sparse([i_indexes;i_indexes+Nmax*(Nmax+2);i_indexes;i_indexes+Nmax*(Nmax+2)],[j_indexes;j_indexes;j_indexes+Nmax*(Nmax+2);j_indexes+Nmax*(Nmax+2)],Qv,2*Nmax*(Nmax+2),2*Nmax*(Nmax+2));
RgQ=sparse([i_indexes;i_indexes+Nmax*(Nmax+2);i_indexes;i_indexes+Nmax*(Nmax+2)],[j_indexes;j_indexes;j_indexes+Nmax*(Nmax+2);j_indexes+Nmax*(Nmax+2)],RgQv,2*Nmax*(Nmax+2),2*Nmax*(Nmax+2));

%solve the T-matrix:
T=-(RgQ/Q)';

if verbose
    figure
    imagesc(abs(T));colormap(jet(256));colorbar
    
    figure
    plot(diag(abs(T).^2),'bx');hold on
    if exist('Tmie','var')
        plot(diag(abs(Tmie).^2),'r');hold off
    else
        hold off
    end
    
    S=2*T+speye(size(T));
    figure
    plot(sum(abs(S).^2));ylim([.95,1.05])
end
