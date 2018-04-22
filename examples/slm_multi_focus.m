% Computation of multiple focusses generated using phase/amplitude SLM
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

% Add the toolbox path
addpath('../');

% Close open figures
close all;
ott.change_warnings('off');

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

%% Setup the incident illumination

%X and Y shall always be in normalised units.
x = linspace(-1, 1, 512);
y = x;
[X,Y]=meshgrid(x, y);

% create incident laser beam, pick a Gaussian beam:
incident_mode=ott.utils.lgmode(0,0,sqrt(X.^2+Y.^2)/fillfactor,atan2(Y,X));
polarisation=[1,0];

%% Plot the laser mode shape on the SLM

figure(1)
contourf(X,Y,incident_mode,32)
hold on
plot(cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'w:','linewidth',1)
hold off
title('mode shape on SLM and target NA zone (in white circle)')
xlabel('X [a.u.]'); ylabel('Y [a.u.]');
axis image

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
xlabel('X [pixels]'); ylabel('Y [pixels]');
axis image
subplot(1,2,2)
imagesc(phase_function)
title('phase of grating')
xlabel('X [pixels]'); ylabel('Y [pixels]');
axis image

%% Compute beam that comes off SLM

% Calculate Nmax: should be beam waist+max translation
Nmax=ott.utils.ka2nmax(2*pi*nMedium*(.5+max(abs(g(:)))));

% Apply the grating to the incident beam
output_mode=G.*incident_mode;

% Calculate the beam shape coefficients of focussed beam
beam = ott.BscPmParaxial(NA, output_mode, 'index_medium', nMedium, ...
    'polarisation', polarisation, 'Nmax', Nmax);

%% create image of the resulting beams along two axes:

[X1,Y1,Z1]=meshgrid(linspace(-3,3)/nMedium,linspace(-3,3)/nMedium,0);
[X2,Z2,Y2]=meshgrid(linspace(-3,3)/nMedium,linspace(-3,3)/nMedium,0);

E = beam.emFieldXyz([[X1(:),Y1(:),Z1(:)];[X2(:),Y2(:),Z2(:)]].');
E2=sum(abs(E).^2,1);

%% plot beams
%NOTE: Due to a legacy choice the plotted field components flip with
%electromagnetic_field_xyz.m. It should be possible to hack your version of
%spharm or electromagnetic_field_xyz.m to fix this field flipping issue.

figure(3)
subplot(1,2,1)
contourf(X1,Y1,reshape(E2(1:end/2),size(X1)),32,'edgecolor','none');
hold on
plot(xlim,mean(ylim)*[1,1],'w:','linewidth',2)
hold off
title('XY fields around focal plane');
xlabel('x [\lambda_m]')
ylabel('y [\lambda_m]')
axis image
subplot(1,2,2)
contourf(X2,Z2,reshape(E2(end/2+1:end),size(X1)),32,'edgecolor','none');
axis image
xlabel('x [\lambda_m]')
ylabel('z [\lambda_m]')
title('XZ fields around focal plane');

%% compute some (approximate) forces on a sphere

T = ott.Tmatrix.simple('sphere', 0.5, 'wavelength0', 1.0, ...
    'index_medium', nMedium, 'index_particle', 1.5);

x = [1;0;0]*linspace(-3,3,20)/nMedium;
y = [0;1;0]*linspace(-3,3,20)/nMedium;
z = [0;0;1]*linspace(-3,3,20)/nMedium;

fxyz1 = ott.forcetorque(beam, T, 'position', x);
fxyz2 = ott.forcetorque(beam, T, 'position', y);
fxyz3 = ott.forcetorque(beam, T, 'position', z);

%% Plot figures of forces
% somewhat confusingly as we are moving the beam these plots are at
% locations negative to their spot locations... double negative fixes this.
figure(4)
subplot(1, 3, 1);
plot(x(1, :),fxyz1.')
hold on
plot(xlim,0*[1,1],'k','linewidth',1)
hold off
grid on
xlabel('x [\lambda_m]')
ylabel('force [Q]')
legend('fx','fy','fz')
title('X-translation');

subplot(1, 3, 2);
plot(y(2, :),fxyz2.')
hold on
plot(xlim,0*[1,1],'k','linewidth',1)
hold off
grid on
xlabel('y [\lambda_m]')
ylabel('force [Q]')
legend('fx','fy','fz')
title('Y-translation');

subplot(1, 3, 3);
plot(z(3, :),fxyz3.')
hold on
plot(xlim,0*[1,1],'k','linewidth',1)
hold off
grid on
xlabel('z [\lambda_m]')
ylabel('force [Q]')
legend('fx','fy','fz')
title('Z-translation');
