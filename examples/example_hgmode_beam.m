% Example of calculation of a HG beam using LG modes in the toolbox.
%
% How long should this take?
% The paraxial part is nearly instantaneous. The part where non-paraxial
% fields are calculated should take quite some time due to grid spacing.
%
% Is there any difficulty?
% Our truncation for LG beams is for the L modes only. We therefore
% recommend underfilling the aperture greatly. It is also vital that all
% the beams have the same rayleigh range so the same w0 must be used for
% all modes.
% 
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Make warnings less obtrusive
ott_warning('once');
change_warnings('off');

%% set up
%first off. which mode do we want?
% we want hg_(3,2). what's the order?
m=3;
n=2;
order=m+n; %order for HG beam.
%order=2*p+l; %order for LG beam.

convergence_angle = 15; % we want a small convergence angle for now.
polarisation = [1,1i];

%% get mode weights
%generate paraxial conversion matrix:
[modeweights,lg_modes,hg_modes]=paraxial_transformation_matrix(order,0,1,0);

%row is the resulting mode
[row]=find(hg_modes(:,1)==m,1);

%let's draw paraxial modes to check!
% note if you do a profile in Z you will need to re-write the code from...
%% here
x=linspace(-pi,pi,128)*45/convergence_angle;
y=x;
z=0;

[X,Y,Z]=meshgrid(x,y,z);

R=sqrt(X.^2+Y.^2+Z.^2);
PHI=atan2(Y,X);

%generate the scalar mode:
U=zeros(size(X));

for ii=1:order+1
    U=U+modeweights(row,ii)*lgmode(lg_modes(ii,1),lg_modes(ii,2),R,PHI,0*R,convergence_angle);
end

UHG=hgmode(m,n,X,Y,0*X,convergence_angle);

h=figure(1);
set(h,'position',[80,80,750,750])
subplot(2,2,1)
imagesc(abs(U).^2);axis equal 
title(['|HG_{' num2str(hg_modes(row,1)) ',' num2str(hg_modes(row,2)) '}|^2' ])
subplot(2,2,2)
imagesc(angle(U));caxis([-pi,pi]);axis equal 
title(['arg(HG_{' num2str(hg_modes(row,1)) ',' num2str(hg_modes(row,2)) '})' ])
subplot(2,2,3)
imagesc(abs(UHG).^2);axis equal 
title(['|HG_{' num2str(m) ',' num2str(n) '}|^2' ])
subplot(2,2,4)
imagesc(angle(UHG));caxis([-pi,pi]);axis equal 
title(['arg(HG_{' num2str(m) ',' num2str(n) '})' ])

% to here
%% If that worked we can now work out the hard part...

% we're going to use the "lean and mean" beam code:

% for low convergence angles a large nmax is required.
nmax=50;
a_full=sparse(nmax*(nmax+2),1);
b_full=sparse(nmax*(nmax+2),1);

for ii=1:order+1

    [n,m,a,b]=bsc_pointmatch_farfield(nmax,1,[lg_modes(ii,:) convergence_angle 1 polarisation 90 ],'sintheta');

    nmax=max(nmax,max(n));

    ci=combined_index(n,m);

    a_full(ci)=a_full(ci)+modeweights(row,ii)*a;
    b_full(ci)=b_full(ci)+modeweights(row,ii)*b;

end

a_full(nmax*(nmax+2)+1:end,1)=false;
b_full(nmax*(nmax+2)+1:end,1)=false;

cia=find(a_full);
cib=find(b_full);

ci=union(cia,cib);

[n,m]=combined_index(ci);
%%
Es=electromagnetic_field_xyz(2*pi*[X(:),Y(:),Z(:)],[n;m],[a_full;b_full],[],[]);

E2=sum(abs(Es.Eincident).^2,2);
H2=sum(abs(Es.Hincident).^2,2);
I=reshape(E2+H2,size(X));

%the transversal phase is approimately the sum of the conjugate weights of
%the jones vectors times the Ex and Ey.
Etr=sum(angle(conj(polarisation(1))*Es.Eincident(:,1)+conj(polarisation(2))*Es.Eincident(:,2)),2);
argEtr=reshape(Etr,size(X));

figure(2)
subplot(1,2,1)
imagesc(I);axis equal
title(['|HG_{' num2str(hg_modes(row,1)) ',' num2str(hg_modes(row,2)) '}|^2' ])
subplot(1,2,2)
imagesc(argEtr);caxis([-pi,pi]);axis equal
title(['arg(HG^{tr}_{' num2str(hg_modes(row,1)) ',' num2str(hg_modes(row,2)) '})' ])
