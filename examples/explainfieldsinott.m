%the thing daryl wanted me to do...

verbose=1;

x=linspace(-2,2,200);
y=0;
z=x;

[X,Y,Z]=meshgrid(x,y,z);

xyz=[X(:),Y(:),Z(:)];

n_medium=1;
lambda=1;

n_particle=1.2;
r_particle=1;

exey=[1,0];

Nmax=12;

k=2*pi;
k_medium=2*pi*n_medium/lambda;

%Nmax=ka2nmax(k_medium*max(r_particle)*max(n_particle)/n_medium);
Nmax=round(k_medium*max(r_particle)*max(n_particle)/n_medium+4.05*(k_medium*max(r_particle)*max(n_particle)/n_medium)^(1/3)+2);
torders=2*Nmax+Nmax^2;

%Nmax=sqrt(size(Tlens,1)/2+1)-1;

rv=r_particle/lambda*n_medium;
nv=n_particle/n_medium;

[RV,NV]=meshgrid(rv,nv);

Ttorders=Nmax^2+2*Nmax;

[nn,mm]=combined_index([1:Ttorders]);

%truncate tmatrix

[n1,m1,a0,b0]=bsc_plane(Nmax,1,0,0,exey(1),exey(2));
[a,b,n,m]=make_beam_vector(a0,b0,n1,m1);

fz=electromagnetic_field_xyz(xyz,[n;m],[a;b]);
r=nmax2ka(12)/2/pi;
theta=linspace(0,2*pi,100);

figure(1)
[~,h]=contourf(x,z,reshape(sum(abs(fz.Eincident),2),size(squeeze(X))).',50)
set(h,'edgecolor','none')
ch=colorbar
colormap(hsv(256))
set(ch,'fontsize',16)
hold on
[~,h]=contour(x,z,reshape(sum(abs(fz.Eincident),2),size(squeeze(X))).',1+[-1e-7,1e-7])
set(h,'edgecolor','k','linewidth',1)
plot(r*cos(theta),r*sin(theta),'w','linewidth',2)
hold off
set(gca,'fontsize',16)
axis equal tight
print(1,'amplitudegoodregion.png','-r300','-dpng')

figure(2)
imagesc(x,z,angle(reshape(fz.Eincident(:,1),size(squeeze(X)))).')
hold on
plot(r*cos(theta),r*sin(theta),'w','linewidth',2)
ch=colorbar
set(ch,'fontsize',16)
colormap(jet(256))
hold off
set(gca,'fontsize',16)
axis equal tight
set(gca,'ydir','normal')
print(2,'phasegoodregion.png','-r300','-dpng')

x=linspace(-2,2,400);
y=0;
z=x;

[X,Y,Z]=meshgrid(x,y,z);

xyz=[X(:),Y(:),Z(:)];

%code to calculate the force acting upon a wavelength radius sphere at the origin
n_relative=1.2; %relative refractive index = 1.2
r_particle=1;   %particle radius is 1 wavelength in the medium
beam_angle=77.80;  %beam half angle of 50 degrees
polarization = [1 0];  %circular polarization
beam_offset = [0 0 0]; %sets the destination of the beam focus to O.

Nmax = ka2nmax(2*pi*r_particle)/2; %calculate limit of the expansion

w0 = lg_mode_w0([0 0], beam_angle);

%calculate the beam shape coefficients
[n, m, a0, b0] = bsc_pointmatch_farfield(Nmax, 1, [0 0 w0 1 polarization 90 beam_offset]);
[a, b, n, m] = make_beam_vector(a0, b0, n, m); %pack beam vector to full size

fz1=electromagnetic_field_xyz(xyz,[n;m],[a;b]);

r=nmax2ka(Nmax)/2/pi;
theta=linspace(0,2*pi,100);

figure(3)
[~,h]=contourf(x,z,reshape(sum(abs(fz1.Eincident),2),size(squeeze(X))).',50)
set(h,'edgecolor','none')
ch=colorbar
colormap(jet(256))
set(ch,'fontsize',16)
hold on
[~,h]=contour(x,z,reshape(sum(abs(fz1.Eincident),2),size(squeeze(X))).',1+[-1e-7,1e-7])
set(h,'edgecolor','k','linewidth',1)
plot(r*cos(theta),r*sin(theta),'w','linewidth',2)
plot(r*cos(theta)-1,r*sin(theta),'w--','linewidth',2)
hold off
set(gca,'fontsize',16)
axis equal tight
axloc=axis;
text(axloc(1)-(axloc(2)-axloc(1))/20,axloc(4),'\bfa)','horizontalalignment','right','verticalalignment','bottom','fontsize',16)
print(3,'amplitudegaussiangoodregion.png','-r300','-dpng')
caxlim=caxis;

figure(4)
imagesc(x,z,angle(reshape(fz1.Eincident(:,1),size(squeeze(X)))).')
hold on
plot(r*cos(theta),r*sin(theta),'w','linewidth',2)
plot(r*cos(theta)-1,r*sin(theta),'w--','linewidth',2)
ch=colorbar
set(ch,'fontsize',16)
colormap(jet(256))
hold off
set(gca,'fontsize',16)
axis equal tight
set(gca,'ydir','normal')
axloc=axis;
text(axloc(1)-(axloc(2)-axloc(1))/20,axloc(4),'\bfb)','horizontalalignment','right','verticalalignment','bottom','fontsize',16)
print(4,'phasegaussiangoodregion.png','-r300','-dpng')

[A,B]=translate_z(Nmax,1);
R=z_rotation_matrix(pi/2,0);
D=wigner_rotation_matrix(Nmax,R);

a1=D'*(A*D*a+B*D*b);
b1=D'*(A*D*b+B*D*a);

fz2=electromagnetic_field_xyz(xyz,[n;m],[a1;b1]);

r=nmax2ka(Nmax)/2/pi;
theta=linspace(0,2*pi,100);

figure(5)
[~,h]=contourf(x,z,reshape(sum(abs(fz2.Eincident),2),size(squeeze(X))).',50)
set(h,'edgecolor','none')
ch=colorbar
colormap(jet(256))
caxis(caxlim);
set(ch,'fontsize',16)
hold on
[~,h]=contour(x,z,reshape(sum(abs(fz2.Eincident),2),size(squeeze(X))).',1+[-1e-7,1e-7])
set(h,'edgecolor','k','linewidth',1)
plot(r*cos(theta),r*sin(theta),'w','linewidth',2)
plot(r*cos(theta)+1,r*sin(theta),'w--','linewidth',2)
hold off
set(gca,'fontsize',16)
axis equal tight
axloc=axis;
text(axloc(1)-(axloc(2)-axloc(1))/20,axloc(4),'\bfc)','horizontalalignment','right','verticalalignment','bottom','fontsize',16)
print(5,'amplitudegaussiantranslated.png','-r300','-dpng')

figure(6)
imagesc(x,z,angle(reshape(fz2.Eincident(:,1),size(squeeze(X)))).')
hold on
plot(r*cos(theta),r*sin(theta),'w','linewidth',2)
plot(r*cos(theta)+1,r*sin(theta),'w--','linewidth',2)
ch=colorbar
set(ch,'fontsize',16)
colormap(jet(256))
hold off
set(gca,'fontsize',16)
axis equal tight
set(gca,'ydir','normal')
axloc=axis;
text(axloc(1)-(axloc(2)-axloc(1))/20,axloc(4),'\bfd)','horizontalalignment','right','verticalalignment','bottom','fontsize',16)
print(6,'phasegaussiangoodtranslated.png','-r300','-dpng')