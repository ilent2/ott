% Example of particle moving in a Gaussian beam trap
%
% Includes configuration for sphere, cylinder and cube.
%
% Approximate of Figure 4 in Nieminen et al., Optical tweezers computational
% toolbox, Journal of Optics A 9, S196-S203 (2007)
%
% How long should this take?
% For the example here it took ~140-160 seconds on a Core2 Duo 6600 with
% 6GB of RAM to calculate the t-matrix for a 1 wavelength sided cube. The
% simluation generally takes less time, here is takes 2 seconds.
%
% Note: The t-matrix code for the cube may be rank deficient to a
% high level. This is a limitation of the PM method. Values on the order of
% 10^-4 may be permissible. Appropriatness of results using such a
% t-matrix depends on the experiment you are comparing it to. This
% instability occurs because of a fundamental difficulty associated with
% fitting a square peg in a round hole!
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add the toolbox to the path (assuming we are in ott/examples)
addpath('../');

% Close open figures
close all;

% Make warnings less obtrusive
ott.warning('once');
ott.change_warnings('off');

% Specify refractive indices
n_medium = 1.33;
n_particle = 1.57;

% Specify the wavelength in the freespace [m]
wavelength0 = 1064e-9;
wavelength_medium = wavelength0 / n_medium;

% Specify the size of the cube
width = 1.0*wavelength_medium;

% Select the type of particle (sphere, cylinder or cube)
particle_type = 'cube';

%% Generate beam

tic

beam = ott.BscPmGauss('NA', 1.25, 'polarisation', [ 1i 1 ], ...
    'power', 1.0, 'index_medium', n_medium, 'wavelength0', wavelength0);

disp(['Beam calculation took ' num2str(toc) ' seconds']);

%% Calculate T-matrix

tic

switch particle_type
  case 'sphere'

    % Describe the shape
    shape = ott.shapes.Shape.simple('sphere', width/2);

    % Specify the initial position
    x_initial = [1;1;1]/5;

    % Set a limit on the max time step, needed for stability, empirically found
    dtlim=0.001;
    
    % Optional parameters for drawing shape
    optargs = {};

  case 'cylinder'

    % Describe the shape
    shape = ott.shapes.Shape.simple('cylinder', [0.5*width, width]);

    % Specify the initial position
    x_initial = [1;1;1]/2;

    % Set a limit on the max time step, needed for stability, empirically found
    dtlim=0.025;
    
    % Optional parameters for drawing shape
    optargs = {'noendcap', true};

  case 'cube'

    % Describe the shape
    shape = ott.shapes.Shape.simple('cube', width);

    % Specify the initial position
    x_initial = [1;1;1]/2;

    % Set a limit on the max time step, needed for stability, empirically found
    dtlim=0.025;
    
    % Optional parameters for drawing shape
    optargs = {'noendcap', true};

  otherwise
    error('Unsupported particle type');
end

T = ott.Tmatrix.simple(shape, ...
    'index_medium', n_medium, ...
    'index_particle', n_particle, ...
    'wavelength0', wavelength0);

[X, Y, Z] = shape.surf('npoints', 20, optargs{:});

disp(['Calculating T-matrix took ', num2str(toc), ' seconds']);

%% Dynamics Simulation
% This dynamics simulation simulates the trap in some hypothetical
% substance. It does not represent any physical system except that it
% displays the same kind of dynamics as a trapped particle capable of
% rotation.

% Time the simulation
tic;

numt=75;            % number of time steps
x=zeros(3,numt);    % array of particle positions
tvec=zeros(numt,1); % array of times where positions are stored

x(:,1)=x_initial;   % initial position for particle

% Set hypothetical drag tensors for the particle
translation_drag_tensor=eye(3)/200;
rotation_drag_tensor=eye(3)/500;

%save rotation matricies for drawing.
Rtotal=zeros(numt*3,3);

%set initial orientation [degrees].
theta = 0;
phi = 0;
Rw = ott.utils.rotz(theta)*ott.utils.roty(phi);
Rtotal(1:3,:) = Rw;

for ii=2:numt

    % Calculate the force and the torque on the particle
    [ft,tt] = ott.forcetorque(beam, T, ...
        'position', x(:, ii-1) * wavelength_medium, 'rotation', Rw);

    % Dynamic time-stepping asymptotic with dtlim. We assume that no
    % multiplier is needed on the rotation to correct the error. There is
    % almost certainly a physically motivated choice for dt, but this method
    % seems to work well in a number of situations. Scaling factors can be
    % added if you care to, let us know what you find out!
    dt=dtlim/(1+(sqrt(sum((inv(rotation_drag_tensor)*tt).^2))));

    % Calculate new x
    x(:,ii)=x(:,ii-1)+inv(translation_drag_tensor)*ft*dt;

    % New rotation matrix using the deviation
    Rw=ott.utils.rotation_matrix(inv(rotation_drag_tensor)*tt*dt)*Rw;
    Rtotal(3*(ii-1)+[1:3],:)=Rw;

    % Store the time
    tvec(ii)=tvec(ii-1)+dt;
end

% Finish timing the calculation
disp(['Simulating dynamics took ' num2str(toc) ' seconds']);

% Create a plot of the trajectory
figure(2)
plot(tvec,x.')
legend('x','y','z')
title('Particle trajectory in principal axes.')
xlabel('Time [simulation units]');
ylabel('Position [\lambda_m]');

%% Plotting code of rotating and translating sphere

x = x * wavelength_medium;

if exist('movieframe', 'var')
  clear movieframe;
end

for ii=1:numt
  figure(1)

  XYZt=Rtotal(3*(ii-1)+[1:3],:)*[X(:).';Y(:).';Z(:).'];
  Xt = reshape(XYZt(1, :) - x(1,ii), size(X));
  Yt = reshape(XYZt(2, :) - x(2,ii), size(Y));
  Zt = reshape(XYZt(3, :) - x(3,ii), size(Z));

  h = surf(Xt, Yt, Zt);

  hold on

  plot3([-1,1],[0,0],[0,0],'b');
  plot3([0,0],[-1,1],[0,0],'g');
  plot3([0,0],[0,0],[-1,1],'r');

  hold off
  xlabel('x');
  ylabel('y');
  zlabel('z');
  grid on
  axis equal
  axis([-1,1,-1,1,-1,1]*width)
  title('Particle trajectory and orientation with time.')
  view(3)
  movieframe(ii)=getframe(1);

  figure(3)

  hold on

  h = surf(Xt, Yt, Zt);
  set(h, 'FaceAlpha', 0);

  if ii==1
    set(h, 'EdgeColor', 'red', 'LineWidth', 2);
  elseif ii==numt
    set(h, 'EdgeColor', 'black', 'LineWidth', 2);
  else
    set(h, 'EdgeColor', 'blue');
  end

  plot3([-1,1],[0,0],[0,0],'b');
  plot3([0,0],[-1,1],[0,0],'g');
  plot3([0,0],[0,0],[-1,1],'r');

  hold off
  xlabel('x');
  ylabel('y');
  zlabel('z');
  grid on
  axis equal
  axis([-1,1,-1,1,-1,1]*width)

  title('Particle trajectory and orientation with time.')
end

%% save movie
if strcmpi(input('Save movie (y/n)? ','s'),'y')
  nm=input('File name (warning overwrites without asking): ','s');
  v = VideoWriter([pwd filesep nm '.avi']);
  open(v);
  writeVideo(v, movieframe);
  close(v);
  clear movieframe;
end

%% Show the figure again

figure(3)
view(2)
