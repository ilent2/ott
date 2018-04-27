% Example of the trapping_landscape code which could be used to produce 
% figure 3 in Nieminen et al., 2007, Journal of Optics. The nuts and bolts 
% of the how is included in landscape() (bellow), this code simply scripts 
% for the parameters used in the article. More analysis of Trapping 
% landscapes and the properties of trapped microspheres appears in Stilgoe 
% et al., 2008, Optics Express. 
%
% The high resolution version takes about ~5000 seconds on a Core2 Duo 6600
% with 6GB RAM. The low resolution version takes about ~1030 seconds on the
% same machine.
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

if exist ("OCTAVE_VERSION", "builtin")
    warning('ott:example_landscape:function',   ...
      ['This code must be modified to run in octave, take the function ', ...
       'defined at the bottom of this script and move it to above where ', ...
       'it is first called.']);
end

% Low res version.
size_range_rad=linspace(1e-2,3.25,100)*1e-6/2; %radius in SI
index_range=linspace(1.34,2.66,50);            %absolute refractive index

% % High res version.
% size_range_rad=linspace(1e-2,3.25,300)*1e-6/2; %radius in SI
% index_range=linspace(1.34,2.66,100);           %absolute refractive index

% if octave plane landscape function immediately below: %



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic
%to see the proceedure for calculating a range of particle sizes and
%refractive indexes open landscape()
n_medium = 1.33;
system_parameters=[1.3 n_medium 1.07e-6 [1 1i] 90 0];
structurelandscape=landscape(index_range,size_range_rad,system_parameters);
toc

%plot the minimum force. This is essentially the figure which appears in
%the article as figure 3.
tempmin=structurelandscape.minforce;
tempmin(tempmin>0)=NaN;
figure(1)
[d,pax]=contourf(size_range_rad*2*1e6,index_range,tempmin,20);
set(pax,'edgecolor','none')
xlabel('particle diameter [{{\mu}{m}}]');
ylabel('refractive index [unitless]');
grid on
cax=colorbar;
ylabel(cax,'Q_z^{min} [unitless]');

%plot the maximum force.
figure(2)
[d,pax]=contourf(size_range_rad*2*1e6,index_range, ...
  structurelandscape.maxforce,20);
set(pax,'edgecolor','none')
xlabel('particle diameter [{{\mu}{m}}]');
ylabel('refractive index [unitless]');
grid on
cax=colorbar;
ylabel(cax,'Q_z^{max} [unitless]');

%plot the z quilibrium position.
figure(3)
[d,pax]=contourf(size_range_rad*2*1e6,index_range,...
  structurelandscape.zequilibrium*n_medium,20);
set(pax,'edgecolor','none')
xlabel('particle diameter [{{\mu}{m}}]');
ylabel('refractive index [unitless]');
grid on
cax=colorbar;
ylabel(cax,'z [k{\lambda}]');

%plot the axial stiffness.
figure(4)
[d,pax]=contourf(size_range_rad*2*1e6,index_range, ...
  structurelandscape.zstiffness/n_medium,20);
set(pax,'edgecolor','none')
xlabel('particle diameter [{{\mu}{m}}]');
ylabel('refractive index [unitless]');
grid on
cax=colorbar;
ylabel(cax,'k_z [k^{-1}\lambda^{-1}]');

%plot the transverse stiffness.
figure(5)
[d,pax]=contourf(size_range_rad*2*1e6,index_range, ...
  structurelandscape.xstiffness/n_medium,20);
set(pax,'edgecolor','none')
xlabel('particle diameter [{{\mu}{m}}]');
ylabel('refractive index [unitless]');
grid on
cax=colorbar;
ylabel(cax,'k_r [k^{-1}\lambda^{-1}]');

function structureOutput=landscape(index_range,size_range_rad,system_parameters)
% landscape.m: Generates a trapping landscape along the beam
%                        propagaing axis.
%
% USAGE:
% [structureOutput]=landscape(index_range,size_range_rad,system_p
% arameters)
%
% VARIABLES:
% index_range, index of materials in medium.
% size_range_rad, radius of particles in medium (um).
% system_parameters, if blank, will calculate for: NA=1.3, no truncation,
% water, wavelength=1 and circular polarization.
% Otherwise, it is a three vector of [NA,medium,lambda,pol,truncation angle (deg)]
% pol is an eigenvector.
% structureOutput=will store the data as a structure. Outputs correctly scaled to RI and
% radius, but not Q or k(Q,\lambda) to power. Re-use the index_range and size_range_rad as plot
% vectors. Z equilibrium position/k will be NaN for any undefined values.

verbose=0;
if nargin<3
    %standard parameters we use... NA 1.3, water, yt fibre, circular pol,90
    %deg truncation (probably should be 77.8).
    system_parameters=[1.3 1.33 1.07e-6 [1 1i] 90 0];
else
    system_parameters= system_parameters(:).';
end

%prep stuff
NA=system_parameters(1);
medium_index=system_parameters(2);
relindex_range=index_range/medium_index;
lambda=system_parameters(3);
relsize_range=size_range_rad/lambda;
k=2*pi*medium_index; %normalized units
nTrans=100;

nParticles=length(relsize_range);
nIndexes=length(relindex_range);
z_equilibrium=zeros(nIndexes,nParticles);
z_stiffness=z_equilibrium;
x_stiffness=z_equilibrium;
maxforcez=z_equilibrium;
minforcez=z_equilibrium;

%% Setup beam

beam = ott.BscPmGauss('NA', NA, 'polarisation', system_parameters(4:5), ...
    'truncation_angle_deg', system_parameters(6), 'power', 1.0, ...
    'k_medium', k, 'index_medium', medium_index);

%% Determine the max Nmax
% This is an optimisation, most of the time you don't need to do this
% since ott.forcetorque does it for you, but the following does
% translations and rotations manually.
maxNmax = ott.utils.ka2nmax(k * max(relsize_range));

%% Run the simulation

for particle=1:nParticles
    disp(['Particle: ' num2str(particle) '/' num2str(nParticles)]);
    if verbose
        disp(['Particle: ' num2str(particle) '/' num2str(nParticles) ', radius=' num2str(relsize_range(particle))]);
    end
    %prepare Nmax and beam vector.
    position=linspace(-(medium_index/NA)^2*(1+relsize_range(particle)*3/4),...
        (medium_index/NA)^2*(1+relsize_range(particle)*3/4),nTrans);
    if relsize_range(particle) > 0.5
        position=linspace(-(medium_index/NA)^2*...
            (0.5+relsize_range(particle)*max(relindex_range)),...
            (medium_index/NA)^2*(0.5+relsize_range(particle)...
            *max(relindex_range)),nTrans);
    end
    %We calculate particle size, then refractive index. Therefore we can
    %calculate all needed translations now.

    % Precompute the beams
    beams = beam.translateZ(position, 'Nmax', maxNmax);

    % Precompute rotation to x-axis
    [~, Dx] = beam.rotateY(pi/2);

    % Precompute the translation on the x-axis (keep A, B invertable)
    dx = 2e-8/beam.k_medium*2*pi;
    [~, A, B] = beam.translateZ(dx/2.0);

    for index=1:nIndexes
        if verbose
        disp(['Index: ' num2str(index) '/' num2str(nIndexes) ', rRI=' num2str(relindex_range(index))]);
        end

        T = ott.Tmatrix.simple('sphere', relsize_range(particle), ...
            'k_medium', k, 'k_particle', k*relindex_range(index));

        % Calculate axial force
        fz = ott.forcetorque(beams, T);

        [maxforcez(index,particle),mai]=max(fz(3, :));
        [minforcez(index,particle),mii]=min(fz(3, :));

        %if minforcez < 1e-8*max(fz)
            if mai<mii && (mii+2) <= nTrans && mai-2 >= 1
                %Splines have good convergence properties within a defined
                %region. This is fast.
                splinex=linspace(position(mai-2),position(mii+2),300);
                splinef=spline(position(mai-2:mii+2),fz(3, mai-2:mii+2));
                spliney=ppval(splinef,splinex);

                maxforcez(index,particle)=max(spliney);
                minforcez(index,particle)=min(spliney);
            else
                %warning('Maximum or minimum doens''t exist!')
                spliney=[];
                splinex=[];
            end

            %Calculate axial equilibrium accurately. I wouldn't trust k for
            %extreme marginal traps anyway.
            if mai<mii && (mii+2) <= nTrans && mai-2 >= 1
                if numel(spliney)>0
                    zero_spline=find(spliney<0,1,'first');
                    if numel(zero_spline)>0
                        %first guess
                        fs=polyfit(splinex(zero_spline-3:zero_spline+2),spliney(zero_spline-3:zero_spline+2),3);

                        rts=roots(fs);
                        rts=rts(imag(rts)==0);

                        [~,zmin]=min(abs(rts-splinex(zero_spline)));

                        z_equilibrium(index,particle)=rts(zmin); %also pretty good.
                        z_stiffness(index,particle)=-polyval([3*fs(1),2*fs(2),fs(3)],rts(zmin)); %not perfect but pretty good.

%                         figure(4)
%                         plot(z_equilibrium(index,particle),0,'r.')
%                         plot([z_equilibrium(index,particle)-position(2)+position(1),z_equilibrium(index,particle)+position(2)-position(1)],-[z_stiffness(index,particle)*(z_equilibrium(index,particle)-position(2)+position(1)),z_stiffness(index,particle)*(z_equilibrium(index,particle)+position(2)-position(1))]+mean([z_stiffness(index,particle)*(z_equilibrium(index,particle)-position(2)+position(1)),z_stiffness(index,particle)*(z_equilibrium(index,particle)+position(2)-position(1))]),'r')
%                         plot([position(1),position(end)],[0,0])
%                         hold off
%                         pause(1)

                        %calc transverse (keep Nmax the same)
                        tbeam = beam.translateZ(z_equilibrium(index,particle));
                        tbeam = tbeam.rotate('wigner', Dx);

                        % Translate along the x direction and calc x force
                        xbeam = tbeam.translate(A, B);
                        [~,~,fx1,~,~,~] = ott.forcetorque(xbeam, T);

                        % Translate along the x direction and calc x force
                        xbeam = tbeam.translate(A', B');
                        [~,~,fx0,~,~,~] = ott.forcetorque(xbeam, T);

                        x_stiffness(index,particle)=(fx1-fx0)/dx;
                    else
                        z_equilibrium(index,particle)=NaN; %also pretty good.
                        z_stiffness(index,particle)=NaN; %not perfect but pretty good.
                        x_stiffness(index,particle)=NaN;
                    end
                end
            else
                z_equilibrium(index,particle)=NaN;
                z_stiffness(index,particle)=NaN;
                x_stiffness(index,particle)=NaN;
            end

    end

end
disp('Output is in calculation units which should be easily convertible into SI.')
structureOutput.maxforce=maxforcez;
structureOutput.minforce=minforcez;
structureOutput.zequilibrium=z_equilibrium;
structureOutput.zstiffness=z_stiffness;
structureOutput.xstiffness=x_stiffness;

end
