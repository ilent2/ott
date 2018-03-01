% Example of calculation of a cube in a Gaussian beam trap
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
% Note for Octave: It is unlikely that the plotting functions at the end
% will work. Just comment out the plotting at the end and replace with the
% Octave equivalent (or nearest) to get some movies! We also ask that if
% you do this that you could send the code to the authors so that it may be
% implemented. Thank you!
%
% PACKAGE INFO

%Tries to clear T matrix if 1. It is good to make 0 if you want repeated
%calculations of the same particle.
clearT=1;

% Specify refractive indices
n_medium = 1.33;
n_particle = 1.57;
n_relative = n_particle/n_medium;

% If you want to give all measurements in wavelengths in the surrounding
% medium, then:
wavelength = 1;
% wavelength = wavelength0 / n_medium;
% else you can give it in any units you want. Only k times lengths matters
k = 2*pi/wavelength;

radius = 1; %this is now side length
Nmax = ka2nmax(k*radius); %nmax for a centered cube with side length radius...
Nmax_medium=Nmax; %these are the same using this method.
Nmax_particle=ka2nmax(k*radius*n_relative); %this is the internal refractive index in this system

k_particle = 2*pi/wavelength*n_relative;

diam_microns = radius * 1.064 * 2 / n_medium;

% Specify the beam width. We can either start with the numerical
% aperture (NA) or the beam convergence angle. Either way, we convert
% to the equivalent paraxial beam waist, which is the w0 we put into the
% paraxial beam to obtain the desired (non-paraxial) far field.
% For a Gaussian beam: w0 = 2/(k*tan(theta))
NA = 1.25;
beam_angle = asin(NA/n_medium)*180/pi;
w0 = lg_mode_w0( [ 0 0 ], beam_angle );

% Polarisation. [ 1 0 ] is plane-polarised along the x-axis, [ 0 1 ] is
% y-polarised, and [ 1 -i ] and [ 1 i ] are circularly polarised.
polarisation = [ i 1 ];

% Location of the focal point relative to the particle. These are the
% [ x y z ] coordinates.
beam_offset = [ 0 0 0];

%Makes beam.
[n,m,a0,b0] = bsc_pointmatch_farfield(Nmax,1,[ 0 0 w0 1 polarisation 90 beam_offset ]);
[a,b,n,m] = make_beam_vector(a0,b0,n,m);

%Finds the root power (note this is different than dividing at the level of
%the forcetorque calculation where we work out squared quantities).
pwr=sqrt(sum(abs(a).^2+abs(b).^2));

%Normalizes force/torque to 1.
a=a/pwr;
b=b/pwr;

%% Insert tmatrix here %%
tic
disp('Calculating T-matrix for cube...')
if clearT | ~exist('T','var')
    T=tmatrix_pm_cube(Nmax,Nmax_medium,Nmax_particle,k,k_particle,radius);
end
% T=tmatrix_mie(Nmax,k,k_particle,radius); %test dynamic simulation with sphere
disp(['It took: ', num2str(toc), ' seconds!'])

%% Dynamics Simulation %%
% This dynamics simulation simulates the trap in some hypothetical
% substance. It does not represent any physical system except that it
% displays the same kind of dynamics as a trapped particle capable of
% rotation.

%Initial conditions:
numt=75; %number of time steps
x=zeros(3,numt); %zeros the length of simulation.
x(:,1)=[1;1;1]/2;
dtlim=0.025; %dt limit (needed for stability near round off). Empirically found.
tvec=zeros(numt,1);

%This is NYI as it involves a model of fluid outside the purview of the
%toolbox. Some future student or user contribution?!?!?
%
%However, in the absense of an interesting shape it does the
%trick of emulating the tranlational and rotational drag of a rotationally
%symetric object. Given the initial values, the normalization factors were
%chosen empirically to give a more interesting simluation.
translation_drag_tensor=eye(3)/200;
rotation_drag_tensor=eye(3)/500; % we assume big rotations. This corresponds to low resistance.

%calculate the x and y-axis rotations for force calculation.
Rx = z_rotation_matrix(pi/2,0);
Dx = wigner_rotation_matrix(Nmax,Rx);

Ry = z_rotation_matrix(pi/2,pi/2);
Dy = wigner_rotation_matrix(Nmax,Ry);

%save rotation matricies for drawing.
Rtotal=zeros(numt*3,3);

%set initial orientation.
Rw=z_rotation_matrix(0*pi/4,0*pi/4);
Rtotal([1:3],:)=Rw;

%set temp force and torque storage
ft=zeros(3,1);
tt=ft;

for ii=2:numt
    rtp=xyz2rtp(x(:,ii-1).').';           %change the coordinate of the old particle position from xyz to spherical polars.
    R = z_rotation_matrix(rtp(2),rtp(3)); %calculates an appropriate axis rotation onto z'.
    D = wigner_rotation_matrix(Nmax,R);   %passive (coordinate) rotation pointing the new particle position along the z axis.
    D2 = wigner_rotation_matrix(Nmax,Rw); %passive rotation of the beam to the particle orientation coordinates.
    
    [A,B] = translate_z(Nmax,rtp(1));     %translation coefficients in the z' direction.
    
    % Rotation of the beam and translation to the new particle location
    % using (D,A,B). D' returns a, b to the original beam coordinates. D2 
    % rotates the beam into particle coordinates.
    a2 = D2*(  D'*A * D*a +  D'*B * D*b ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
    b2 = D2*(  D'*A * D*b +  D'*B * D*a ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
    
    % calculate the scattered light in the particle coordinates.
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    %Return to from particle coordinates original beam coordinates. Note: we
    %can apply Rw to the force and torque vectors instead as it is an
    %active rotation (particle rotation) whereas D2 is a passive rotation
    %(coordinate rotation).
    ad=D2'*a2;
    bd=D2'*b2;
    pd=D2'*p;
    qd=D2'*q;
    
%     [ft(3),tt(3)] = force_z(n,m,ad,bd,pd,qd);             %z-force calculation.
%     [ft(1),tt(1)] = force_z(n,m,Dx*ad,Dx*bd,Dx*pd,Dx*qd); %Dx makes the z-force calculation the x-force calculation.
%     [ft(2),tt(2)] = force_z(n,m,Dy*ad,Dy*bd,Dy*pd,Dy*qd); %Dy makes the z-force calculation the y-force calculation.
    
    %This method is more streamlined than the above.
    [ft,tt] = force_torque_farsund(n,m,ad,bd,pd,qd);
    
    %Dynamic time-stepping asymptotic with dtlim. We assume that no
    %multiplier is needed on the rotation to correct the error. There is
    %almost certainly a physically motivated choice for dt, but this method
    %seems to work well in a number of situations. Scaling factors can be
    %added if you care to, let us know what you find out!
    dt=dtlim/(1+(sqrt(sum((-inv(rotation_drag_tensor)*tt).^2))));
    
    %Calculate new x.
    x(:,ii)=x(:,ii-1)+inv(translation_drag_tensor)*ft*dt;
    
    % %These are examples of using Rw instead:
    %x(:,ii)=x(:,ii-1)+inv(Rw'*translation_tensor)*ft*dt;
    %x(:,ii)=x(:,ii-1)+inv(translation_tensor)*Rw*ft*dt;
    % %You may find using the above quicker than inverting D2. This
    % %involves going up several lines and changing ad, bd, pd, qd back to
    % %a2, b2, p, q.
    
    %New rotation matrix using the deviation.
    Rw=rotation_matrix(-inv(rotation_drag_tensor)*tt*dt)*Rw; %This calculates the appropriate rotation using the euler-rodriguez formula.
    Rtotal(3*(ii-1)+[1:3],:)=Rw;
    
    tvec(ii)=tvec(ii-1)+dt; %update time
end

figure(2)
plot(tvec,x.')
legend('x','y','z')
title('Particle trajectory in principal axes.')

%% Plotting code of rotating and translating sphere %%

% patch of cube
verticies=[-1,-1,-1;1,-1,-1;1,1,-1;-1,1,-1;-1,-1,1;1,-1,1;1,1,1;-1,1,1]*radius/2;
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
        XYZt=Rtotal(3*(ii-1)+[1:3],:)*[verticies(jj,1);verticies(jj,2);verticies(jj,3)];
        Xt(jj)=XYZt(1);
        Yt(jj)=XYZt(2);
        Zt(jj)=XYZt(3);
    end
    
    if ii==1
        patch('Vertices',[Xt-x(1,ii),Yt-x(2,ii),Zt-x(3,ii)],'Faces',faces,'FaceVertexCData',hsv(6),'FaceColor','flat','edgecolor','red');
    elseif ii==numt
        patch('Vertices',[Xt-x(1,ii),Yt-x(2,ii),Zt-x(3,ii)],'Faces',faces,'FaceVertexCData',hsv(6),'FaceColor','flat','edgecolor','black');
    else
        patch('Vertices',[Xt-x(1,ii),Yt-x(2,ii),Zt-x(3,ii)],'Faces',faces,'FaceVertexCData',hsv(6),'FaceColor','flat','edgecolor','blue');
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
    axis([-1,1,-1,1,-1,1])
    title('Particle trajectory and orientation with time.')
    view(3)
    movieframe(ii)=getframe(1);
    
    figure(3)
    if ii==1
        patch('Vertices',[Xt-x(1,ii),Yt-x(2,ii),Zt-x(3,ii)],'Faces',faces,'Facealpha',0,'edgecolor','red','linewidth',2);
    elseif ii==numt
        patch('Vertices',[Xt-x(1,ii),Yt-x(2,ii),Zt-x(3,ii)],'Faces',faces,'Facealpha',0,'edgecolor','black','linewidth',2);
    else
        patch('Vertices',[Xt-x(1,ii),Yt-x(2,ii),Zt-x(3,ii)],'Faces',faces,'Facealpha',0,'edgecolor','blue');
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

%% save movie %%
if strcmp(lower(input('Save movie (y/n)? ','s')),'y')
    nm=input('File name (warning overwrites without asking): ','s');
    upath=userpath;
    movie2avi(movieframe,[upath(1:end-1),'\',nm,'.avi']);
% else
%     try; close 1;end;
%     figure(1)
%     movie(1,movieframe,1,30);
end

figure(3)
view(2)