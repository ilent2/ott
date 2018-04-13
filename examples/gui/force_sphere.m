function [x,fx,y,fy,z,fz]=force_sphere(Nmax,particle_params,beam_params,scanrange,handles)
%% Initial parameters
% For testing
% Nmax=12;
% r=1.2;
% n_particle=1.51;
% n_medium=1.33;

% wavelength=1;
% polarisation=[1 0];
% NA =1.2;
% beam_angle = asin(NA/n_medium)*180/pi;
% b_mode=[0 0];    % Beam type
% w0 = lg_mode_w0( b_mode, beam_angle );
% beam_offset = [ 0 0 0];
% truncation_angle=90;
% 
% Nx=50;
% Ny=50;
% Nz=50;
% xrange=[-2.5,2.5];
% yrange=[-2.5,-2.5];
% zrange=[-2.5,-2.5];
%% Input parameters

b_type=beam_params.beam_type;
b_mode=beam_params.mode;
k=beam_params.k;
polarisation=beam_params.polarisation;
beam_offset=beam_params.beam_offset;
w0=beam_params.w0;
truncation_angle=beam_params.truncation_angle;

rad=particle_params.Radius;
n_particle=particle_params.nparticle;
n_medium=particle_params.nmedium;
n_relative = n_particle/n_medium;

xrange=scanrange.x;
yrange=scanrange.y;
zrange=scanrange.z;

tlefthdl=handles.timeleft;
tlefttexthdl=handles.timelefttext;
gplothdl=handles.gauge;
xplothdl=handles.xplot;
zplothdl=handles.zplot;

%% Beam and particle coeff.

% [ p l w0 P xcomponent ycomponent truncation_angle xoffset yoffset zoffset ]
[n0,m0,a0,b0] = bsc_pointmatch_farfield(Nmax,1,[ b_mode w0 1 polarisation truncation_angle beam_offset]);
[a,b,n,m]=make_beam_vector(a0,b0,n0,m0,Nmax);

% tmatrix_mie(Nmax,k_medium,k_particle,radius)

T = tmatrix_mie(Nmax,k,k*n_relative,rad);

% z = linspace(zrange(1),zrange(2),Nz);
% x = linspace(xrange(1),xrange(2),Nx);
% y = linspace(yrange(1),yrange(2),Ny);
[x,y,z] = meshgrid(xrange,yrange,zrange);

% 
% fx = zeros(size(x,2),size(y,2),size(z,2));
% fy = zeros(size(x,2),size(y,2),size(z,2));
% fz = zeros(size(x,2),size(y,2),size(z,2));

fx = zeros(1,1,1);
fy = zeros(1,1,1);
fz = zeros(1,1,1);

%root power for nomalization to a and b individually.
pwr = sqrt(sum( abs(a).^2 + abs(b).^2 ));
%normalize total momentum of wave sum to 1. Not good for SI EM field.
a=a/pwr;
b=b/pwr;

cmp=[linspace(1,0,101)',linspace(1,1,101)',linspace(1,0,101)'];
% cmp=flip(hot(200),1);
gplothdl.Value=0;
for i=1:size(yrange,2)
    for j=1:size(xrange,2)
        tic;
        for l=1:size(zrange,2)
            [rt,theta,phi]=xyz2rtp(x(i,j,l),y(i,j,l),z(i,j,l));
           
            [a2,b2,p2,q2]=coef3DTranslation(a,b,T,rt,theta,phi,Nmax);
            [fxi,fyi,fzi,txi,tyi,tzi]=force_torque_farsund(n,m,a2,b2,p2,q2);
            fx(i,j,l)=fxi;
            fy(i,j,l)=fyi;
            fz(i,j,l)=fzi;
        end
        gplothdl.Value=(j+(i-1)*size(yrange,2))*100/(size(yrange,2)*size(xrange,2));
        t=toc;
        tlefttexthdl.Value=datestr(t*(size(yrange,2)*size(xrange,2)-(j+(i-1)*size(yrange,2)))/(24*60*60), 'HH:MM:SS.FFF');
        gplothdl.BackgroundColor=cmp(ceil((j+(i-1)*size(yrange,2))*100/(size(yrange,2)*size(xrange,2))),:);
    end
end
% figure,mesh(fx(:,:,5),fz(:,:,5))
% [xplot,yplot,zplot]=meshgrid(x,y,z);
% figure,quiver3(xplot,yplot,zplot,fx,fy,fz)
% figure,quiver3(x,y,z,fx,fy,fz),axis equal


% for i=-2.5:0.01:2.5
% %    h=contourslice(x,y,z,sqrt(fx.^2+fy.^2+fz.^2),i,i,i,1);axis equal;
%    h=contourslice(x,y,z,fz,i,i,i,1);axis equal;
% %     h=slice(x,y,z,sqrt(fx.^2+fy.^2+fz.^2),2.5,i,-2.5);
% 
%     drawnow
%     pause(0)
%     
% end
%% Plots

% figure,h=contourslice(x,y,z,fz,-2.5:0.01:2.5,-2.5:0.01:2.5,-2.5:0.01:2.5,1);axis equal;