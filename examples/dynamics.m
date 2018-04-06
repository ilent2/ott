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

% Make warnings less obtrusive
ott_warning('once');
change_warnings('off');

% Specify refractive indices
n_medium = 1.33;
n_particle = 1.57;

% Specify the wavelength in the freespace [m]
wavelength0 = 1064e-9;

% Specify the size of the cube
width = 1.0*wavelength0/n_medium;

% Select the type of particle (sphere, cylinder or cube)
particle_type = 'cube';

%% Generate beam

beam = ott.BscPmGauss('NA', 1.25, 'polarisation', [ i 1 ], 'power', 1.0);

%% Calculate T-matrix

tic
switch particle_type
  case 'sphere'
    T = ott.Tmatrix.simple('sphere', 0.5*width, ...
        'n_medium', n_medium, ...
        'n_particle', n_particle, ...
        'wavelength0', wavelength0);

  case 'cylinder'
    T = ott.Tmatrix.simple('cylinder', [ 0.5*width, width ] ...
        'n_medium', n_medium, ...
        'n_particle', n_particle, ...
        'wavelength0', wavelength0);

  case 'cube'
    T = ott.Tmatrix.simple('cube', width, ...
        'n_medium', n_medium, ...
        'n_particle', n_particle, ...
        'wavelength0', wavelength0);

  otherwise
    error('Unsupported particle type');
end
disp(['Calculating T-matrix took ', num2str(toc), ' seconds!']);

%% Dynamics Simulation
% This dynamics simulation simulates the trap in some hypothetical
% substance. It does not represent any physical system except that it
% displays the same kind of dynamics as a trapped particle capable of
% rotation.

numt=75;            % number of time steps
x=zeros(3,numt);    % array of particle positions
tvec=zeros(numt,1); % array of times where positions are stored

x(:,1)=[1;1;1]/2;   % initial position for particle

% Set a limit on the max time step, needed for stability, empirically found
dtlim=0.025;

% Set hypothetical drag tensors for the particle
translation_drag_tensor=eye(3)/200;
rotation_drag_tensor=eye(3)/500;

%save rotation matricies for drawing.
Rtotal=zeros(numt*3,3);

%set initial orientation.
Rw = rotz(0)*roty(0);
Rtotal([1:3],:) = Rw;

for ii=2:numt

    % Translate the beam to the particle
    tbeam = beam.translateXyz(x(:, ii-1));

    % Rotate the beam to match the particle rotation
    [rbeam, D] = tbeam.rotate(Rw);

    % Scatter the beam and return from the particle frame to the lab frame
    sbeam = D' * (T * rbeam);
    
    % Calculate the force and the torque on the particle
    [ft,tt] = ott.forcetorque(tbeam, sbeam);
    
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

% Create a plot of the trajectory
figure(2)
plot(tvec,x.')
legend('x','y','z')
title('Particle trajectory in principal axes.')

%% Plotting code of rotating and translating sphere

% patch of cube
verticies=[-1,-1,-1;1,-1,-1;1,1,-1;-1,1,-1; ...
    -1,-1,1;1,-1,1;1,1,1;-1,1,1]*radius/2;
faces=[1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
cdata=zeros(size(verticies));

Xt=zeros(length(verticies(:,1)),1);
Yt=Xt;
Zt=Xt;

try; clear movieframe;end;

for ii=1:numt
  figure(1)
  clf
  for jj=1:length(verticies(:,1))
      XYZt=Rtotal(3*(ii-1)+[1:3],:)* ...
          [verticies(jj,1);verticies(jj,2);verticies(jj,3)];
      Xt(jj)=XYZt(1);
      Yt(jj)=XYZt(2);
      Zt(jj)=XYZt(3);
  end
  
  if ii==1
    patch('Vertices',[Xt-x(1,ii),Yt-x(2,ii),Zt-x(3,ii)], ...
        'Faces',faces,'FaceVertexCData',hsv(6), ...
        'FaceColor','flat','edgecolor','red');
  elseif ii==numt
    patch('Vertices',[Xt-x(1,ii),Yt-x(2,ii),Zt-x(3,ii)], ...
        'Faces',faces,'FaceVertexCData',hsv(6), ...
        'FaceColor','flat','edgecolor','black');
  else
    patch('Vertices',[Xt-x(1,ii),Yt-x(2,ii),Zt-x(3,ii)], ...
        'Faces',faces,'FaceVertexCData',hsv(6), ...
        'FaceColor','flat','edgecolor','blue');
  end
  
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
  axis([-1,1,-1,1,-1,1]*radius)
  title('Particle trajectory and orientation with time.')
  view(3)
  movieframe(ii)=getframe(1);
  
  figure(3)
  if ii==1
    patch('Vertices',[Xt-x(1,ii),Yt-x(2,ii),Zt-x(3,ii)], ...
        'Faces',faces,'Facealpha',0,'edgecolor','red','linewidth',2);
  elseif ii==numt
    patch('Vertices',[Xt-x(1,ii),Yt-x(2,ii),Zt-x(3,ii)], ...
        'Faces',faces,'Facealpha',0,'edgecolor','black','linewidth',2);
  else
    patch('Vertices',[Xt-x(1,ii),Yt-x(2,ii),Zt-x(3,ii)], ...
        'Faces',faces,'Facealpha',0,'edgecolor','blue');
  end
  
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
  
  title('Particle trajectory and orientation with time.')
end

%% save movie
if strcmp(lower(input('Save movie (y/n)? ','s')),'y')
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
