%Example of generating a HG beam from LG beam.
%Define a HG mode
m=3;
n=2;

%The order is the sum of the mode indices.
order=m+n;

%For non-paraxial beams radial modes don't appear near the focus...
convergence_angle = 15; % we want a small convergence angle for now.
polarisation = [1,1i];

%generate paraxial conversion matrix:
[modeweights,lglookups,hglookups]=genLG2HG(order);

%we don't have reverse lookup so we're going to find the row of
%the transformation matrix using find.
[m_,n_]=hglookup(order,hglookups);

%row is the resulting mode for the LG to HG conversion.
[row,col]=find(m_==m,1);

%Make a plot grid for the transverse profile.
[X,Y]=meshgrid(linspace(-3,3,128),linspace(-3,3,128));
R=sqrt(X.^2+Y.^2);
PHI=atan2(Y,X);
U=zeros(size(PHI));
for ii=1:order+1
    [p,l]=lglookup(order,lglookups(row,ii));
    U=U+modeweights(row,ii)*lgmode(p,l,R,PHI);
end

%Generate the hgmode for comparison:
UHG=hgmode(m,n,X,Y);

%Draw a fancy plot of the intensity and phase.
h=figure(1)
set(h,'position',[80,80,750,750])
subplot(2,2,1)
imagesc(abs(U).^2);axis equal 
title(['LG derived |HG_{' num2str(m_(row,1)) ',' num2str(n_(row,1)) '}|^2' ],'fontsize',16)
subplot(2,2,2)
imagesc(angle(U));caxis([-pi,pi]);axis equal 
title(['LG derived arg(HG_{' num2str(m_(row,1)) ',' num2str(n_(row,1)) '})' ],'fontsize',16)
subplot(2,2,3)
imagesc(abs(UHG).^2);axis equal 
title(['|HG_{' num2str(m) ',' num2str(n) '}|^2' ],'fontsize',16)
subplot(2,2,4)
imagesc(angle(UHG));caxis([-pi,pi]);axis equal 
title(['arg(HG_{' num2str(m) ',' num2str(n) '})' ],'fontsize',16)

% as bsc_lgmode_farfield.m will expand to nmax of 200
% let's set a's and b's to the full size.
a_full=sparse(200*202,1);
b_full=sparse(200*202,1);

for ii=1:order+1
    
    [p,l]=lglookup(order,lglookups(row,ii));

    [n,m,a,b]=bsc_lgmode_farfield(convergence_angle,[p,l],'fixed',polarisation);
    
    cit=combined_index(n,m);
    ci=unique([ci;cit]);

    a_full(cit)=a_full(cit)+modeweights(row,ii)*a;
    b_full(cit)=b_full(cit)+modeweights(row,ii)*b;

end
%Let's now find the non-zero n,m elements...
[n,m]=combined_index(ci);

%Truncate to nmax
a_full(max(n)*(max(n)+2)+1:end,1)=false;
b_full(max(n)*(max(n)+2)+1:end,1)=false;

%Generate field
Es=electromagnetic_field_xyz([X(:),Y(:),0*X(:)],[n;m],[a_full;b_full]);
E2=sum(abs(Es.Eincident).^2,2);
I=reshape(E2,size(X));

%the transversal phase is approimately the sum of the conjugate
%weights of the jones vectors times the Ex and Ey.
Etr=sum(conj(polarisation(1))*Es.Eincident(:,1)+conj(polarisation(2))*Es.Eincident(:,2),2);
argEtr=reshape(angle(Etr),size(X));

h=figure(2)
set(h,'position',[100,100,750,375])
subplot(1,2,1)
imagesc(I);axis equal tight
title(['|HG_{' num2str(m_(row,1)) ',' num2str(n_(row,1)) '}|^2' ],'fontsize',16)
subplot(1,2,2)
imagesc(argEtr);caxis([-pi,pi]);axis equal tight
title(['arg(HG^{tr}_{' num2str(m_(row,1)) ',' num2str(n_(row,1)) '})' ],'fontsize',16)