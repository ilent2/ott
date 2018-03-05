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
% PACKAGE INFO

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
[modeweights,lglookups,hglookups]=genLG2HG(order);

%we don't have reverse lookup so we're going to find the appropriate row of
%the transformation matrix using find.
[m_,n_]=hglookup(order,hglookups);

%row is the resulting mode
[row,col]=find(m_==m,1);

%let's draw paraxial modes to check!
% note if you do a profile in Z you will need to re-write the code from...
%% here
x=linspace(-3,3,128);
y=x;
z=0;

[X,Y,Z]=meshgrid(x,y,z);

R=sqrt(X.^2+Y.^2+Z.^2);
PHI=atan2(Y,X);

%generate the scalar mode:
U=zeros(size(X));

for ii=1:order+1
    [p,l]=lglookup(order,lglookups(row,ii));
    U=U+modeweights(row,ii)*lgmode(p,l,R,PHI);
end

UHG=hgmode(m,n,X,Y);

h=figure(1)
set(h,'position',[80,80,750,750])
subplot(2,2,1)
imagesc(abs(U).^2);axis equal 
title(['|HG_{' num2str(m_(row,1)) ',' num2str(n_(row,1)) '}|^2' ])
subplot(2,2,2)
imagesc(angle(U));caxis([-pi,pi]);axis equal 
title(['arg(HG_{' num2str(m_(row,1)) ',' num2str(n_(row,1)) '})' ])
subplot(2,2,3)
imagesc(abs(UHG).^2);axis equal 
title(['|HG_{' num2str(m) ',' num2str(n) '}|^2' ])
subplot(2,2,4)
imagesc(angle(UHG));caxis([-pi,pi]);axis equal 
title(['arg(HG_{' num2str(m) ',' num2str(n) '})' ])

% to here
%% If that worked we can now work out the hard part...

% we're going to use the "lean and mean" beam code:

% as the bsc_lgmode_farfield will expand to nmax of 200 let's set a's and
% b's to the full size.

a_full=sparse(200*202,1);
b_full=sparse(200*202,1);

nmax=0;
for ii=1:order+1
    
    [p,l]=lglookup(order,lglookups(row,ii));

    [n,m,a,b]=bsc_lgmode_farfield(convergence_angle,[p,l],'fixed',polarisation);

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

Es=electromagnetic_field_xyz([X(:),Y(:),Z(:)],[n;m],[a_full;b_full]);

E2=sum(abs(Es.Eincident).^2,2);
I=reshape(E2,size(X));

%the transversal phase is approimately the sum of the conjugate weights of
%the jones vectors times the Ex and Ey.
Etr=sum(angle(conj(polarisation(1))*Es.Eincident(:,1)+conj(polarisation(2))*Es.Eincident(:,2)),2);
argEtr=reshape(Etr,size(X));

figure(2)
subplot(1,2,1)
imagesc(I);axis equal
title(['|HG_{' num2str(m_(row,1)) ',' num2str(n_(row,1)) '}|^2' ])
subplot(1,2,2)
imagesc(argEtr);caxis([-pi,pi]);axis equal
title(['arg(HG^{tr}_{' num2str(m_(row,1)) ',' num2str(n_(row,1)) '})' ])