% example_slm.m computation for SLM.
% 
% Computes the BSCs for a three beam optical tweezers generated with phase 
% only SLM. Instead of using a gratings and lenses algorithm, it translates
% using angular phase shift approximations computed from the expression:
% exp(1i*dot(k,r));
%
% It should generate the correct beam shift for all of the beams generated.
% The near fields are plotted at the end.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

%define microscope target NA:
NA=-1.33*sin(64.46*pi/180); %NA in water... (from bottom)
nMedium=1.33;

% Due to issues with numerical stability, compute G&L around first order:
% Create "three beam" grating function:
g=[2,0,0;0,2,0;0,0,2]; %
% generate "l" modes:
l=[0;0;0]; %

% set up SLM bounds. To underfill change fillfactor to a number less than 1.
fillfactor=1;

% Is SLM in phase and or amplitude modes?
phase_mode=true;
amplitude_mode=false;

%X and Y shall always be in normalised units.
[X,Y]=meshgrid(linspace(-1,1,512),linspace(-1,1,512));

% create incident laser beam, pick a Gaussian beam:
incident_mode=lgmode(0,0,sqrt(X.^2+Y.^2)/fillfactor,atan2(Y,X));
% incident_mode=hgmode(0,0,X,Y);
polarisation=[1,0];

%% Plot the laser mode shape on the SLM

figure(1)
contourf(X,Y,incident_mode,32)
hold on
plot(cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'w:','linewidth',1)
hold off
title('mode shape on SLM and target NA zone (in white circle)')
axis equal

%% Create a pattern for the SLM
%NOTE: This is not actually a paraxial system. We are assuming that
%this field is proportional to the filling r and not the angle theta
%create actual grating first and then compute SLM operation mode:
G=zeros(size(X));
%This is a rescaling step for the purposes of elegance... guaranteed:
XG=X*NA/nMedium;
YG=Y*NA/nMedium;
PHI=atan2(YG,XG)/2/pi; %phi normalised from [-pi, pi] to [-.5,.5];

% phase target: exp(1i*dot(k,r)). XG/YG are already \propto sin(theta).
% Compute \propto cos(theta) component for z-motion:
costheta=sign(NA)*sqrt(1-XG.^2-YG.^2);
costheta(logical(imag(costheta)))=0;
for ii=1:size(g,1)
    G=G+exp(1i*2*pi*(g(ii,1)*XG+g(ii,2)*YG+g(ii,3)*costheta+l(ii)*PHI));
end

phase_function=ones(size(X));
amplitude_function=phase_function;

if phase_mode
    phase_function=angle(G);
end
if amplitude_mode
    amplitude_function=abs(G)/max(abs(G(:)));
end

G=amplitude_function.*exp(1i*phase_function);

%% Plot the amplitude and phase of the grating

figure(2)
subplot(1,2,1)
imagesc(amplitude_function)
title('amplitude of grating')
axis equal
subplot(1,2,2)
imagesc(phase_function)
title('phase of grating')
axis equal


%% Compute beam that comes off SLM

% Calculate Nmax: should be beam waist+max translation
Nmax=ka2nmax(2*pi*nMedium*(.5+max(abs(g(:)))));

% Apply the grating to the incident beam
output_mode=G.*incident_mode;

% Calculate the beam shape coefficients of focussed beam
[a,b]=paraxial_to_bsc(NA,output_mode,polarisation,{'nmax',Nmax});

%% create image of the resulting beams along two axes:

[X1,Y1,Z1]=meshgrid(2*pi*linspace(-3,3),2*pi*linspace(-3,3),0);
[X2,Z2,Y2]=meshgrid(2*pi*linspace(-3,3),2*pi*linspace(-3,3),0);

[n,m]=combined_index([1:length(a)].');
Es=electromagnetic_field_xyz([[X1(:),Y1(:),Z1(:)];[X2(:),Y2(:),Z2(:)]],[n;m],[a;b]);

E2=sum(abs(Es.Eincident).^2,2);

%% plot beams
%NOTE: Due to a legacy choice the plotted field components flip with
%electromagnetic_field_xyz.m. It should be possible to hack your version of
%spharm or electromagnetic_field_xyz.m to fix this field flipping issue.

figure(3)
subplot(1,2,1)
contourf(X1/2/pi,Y1/2/pi,reshape(E2(1:end/2),size(X1)),32,'edgecolor','none');
hold on
plot(xlim,mean(ylim)*[1,1],'w:','linewidth',2)
hold off
title('XY fields around focal plane');
xlabel('x')
ylabel('y')
axis equal
subplot(1,2,2)
contourf(X2/2/pi,Z2/2/pi,reshape(E2(end/2+1:end),size(X1)),32,'edgecolor','none');
axis equal
xlabel('x')
ylabel('z')
title('XZ fields around focal plane');

%% compute some (approximate) forces on a sphere

T=tmatrix_mie(Nmax,2*pi*nMedium,2*pi*1.5,.5);

x=linspace(-3,3,20);
y=linspace(-3,3,20);
z=linspace(-3,3,20);

fxyz1=zeros(3,length(x));
fxyz2=zeros(3,length(y));
fxyz3=zeros(3,length(z));

[rtp1]=xyz2rtp(1*x,0*x,0*x);
[rtp2]=xyz2rtp(0*y,1*y,0*y);
[rtp3]=xyz2rtp(0*z,0*z,1*z);

%calculate the force along x,y,z independently
for ii = 1:length(z)
    
    [A,B] = translate_z(Nmax,rtp3(ii,1)); % use same translation

    W=wigner_rotation_matrix(Nmax,rotation_matrix([0,1,0],(1-sign(z(ii)))*pi/2));

    a2 = W'*( A*W*a + B*W*b );
    b2 = W'*( A*W*b + B*W*a );
    
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    fxyz3(:,ii) = forcetorque(n,m,a2,b2,p,q);

    W=wigner_rotation_matrix(Nmax,rotation_matrix([-sin(rtp1(ii,3)),cos(rtp1(ii,3)),0],rtp1(ii,2)));

    a2 = W'*( A*W*a + B*W*b );
    b2 = W'*( A*W*b + B*W*a );
    
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    fxyz1(:,ii) = forcetorque(n,m,a2,b2,p,q);

    W=wigner_rotation_matrix(Nmax,rotation_matrix([-sin(rtp2(ii,3)),cos(rtp2(ii,3)),0],rtp2(ii,2)));

    a2 = W'*( A*W*a + B*W*b );
    b2 = W'*( A*W*b + B*W*a );
    
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    fxyz2(:,ii) = forcetorque(n,m,a2,b2,p,q);
    
end
%% Plot figures of forces
% somewhat confusingly as we are moving the beam these plots are at
% locations negative to their spot locations... double negative fixes this.
figure(4)
subplot(1, 3, 1);
plot(x,fxyz1.')
hold on
plot(xlim,0*[1,1],'k','linewidth',1)
hold off
grid on
xlabel('x displacement')
ylabel('force [Q]')
legend('fx','fy','fz')
title('X-translation');

subplot(1, 3, 2);
plot(y,fxyz2.')
hold on
plot(xlim,0*[1,1],'k','linewidth',1)
hold off
grid on
xlabel('y displacement')
ylabel('force [Q]')
legend('fx','fy','fz')
title('Y-translation');

subplot(1, 3, 3);
plot(z,fxyz3.')
hold on
plot(xlim,0*[1,1],'k','linewidth',1)
hold off
grid on
xlabel('z displacement')
ylabel('force [Q]')
legend('fx','fy','fz')
title('Z-translation');
