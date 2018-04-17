% Example of calculation of restoring forces on a layered sphere.
%
% Approximate of Figure 3a/b in Hu et al., Antireflection coating for 
% improved optical trapping, Journal of Applied Physics 103, 093119 (2008)
%
% How long should this take?
% Using the the parameters required to fairly faithfully reproduce figure 3
% took 800 seconds on a Core2 Duo 6600 w/ 6GB RAM. The version here uses a
% much coarser grid and takes 200 seconds.
%
% Note: This is not fully annotated. example_gaussian.m may be a better
% place to start for those unfamiliar with the toolbox.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Make warnings less obtrusive
ott_warning('once');
change_warnings('off');

% Specify refractive indices
n_medium = 1.34;
n_particle = 1.8;

% Shell refractive index
n_shell = sqrt(n_medium*n_particle);

% If you want to give all measurements in wavelengths in freespace,
% then:
wavelength = 1;
% wavelength = wavelength0 / n_medium;
% else you can give it in any units you want. Only k times lengths matters
k = 2*pi/wavelength;

tic

%radius = linspace(.25,2.5,200)/n_medium; %in freespace units...
radius = linspace(.25,2.5,50)/n_medium; %in freespace units...

Nmax = ka2nmax(k*(radius+1/n_shell/4));

diam_microns = radius * 1.064 * 2;
%
%     if Nmax < 12
%         Nmax = 12;
%     end

% Specify the beam width. We can either start with the numerical
% aperture (NA) or the beam convergence angle. Either way, we convert
% to the equivalent paraxial beam waist, which is the w0 we put into the
% paraxial beam to obtain the desired (non-paraxial) far field.
% For a Gaussian beam: w0 = 2/(k*tan(theta))
NA = 1.02;
beam_angle = asin(NA/n_medium)*180/pi;

% Polarisation. [ 1 0 ] is plane-polarised along the x-axis, [ 0 1 ] is
% y-polarised, and [ 1 -i ] and [ 1 i ] are circularly polarised.
polarisation = [ 1 i ];

[n,m,a0,b0] = bsc_pointmatch_farfield(max(Nmax),1,[ 0 0 beam_angle 1 polarisation 90 ]);
[a,b,n,m] = make_beam_vector(a0,b0,n,m);
%root power for nomalization to a and b individually.
pwr = sqrt(sum( abs(a).^2 + abs(b).^2 ));
%normalize total momentum of wave sum to 1. Not good for SI EM field.
a=a/pwr;
b=b/pwr;

z = linspace(-.5,max(radius)*n_particle,45)*n_medium; %in medium units
r = linspace(0,max(radius)*1.33,25)*n_medium; %in medium units

fz = zeros(size(z));
fr = zeros(size(r));
fz1=fz;
fr1=fr;
axialrestoringforce=zeros([size(radius),1]);
axialrestoringforce1=axialrestoringforce;
transverserestoringforce=zeros([size(radius),1]);
transverserestoringforce1=transverserestoringforce;
for ii=1:length(radius)

    disp(['Calculating radius ' num2str(ii) '/' ...
        num2str(length(radius)) ' (' num2str(radius(ii)) ')']);

    %re-calculate grid because it's more accurate.
    z = linspace(-.5,max([1,radius(ii)])*n_particle,45)*n_medium; %in medium units
    r = linspace(0,max([1,radius(ii)])*1.33,25)*n_medium; %in medium units

    T = tmatrix_mie_layered(max(Nmax),k*n_medium,k*[n_particle,n_shell],([radius(ii),radius(ii)+1/n_shell/4]));
    T1 = tmatrix_mie_layered(max(Nmax),k*n_medium,k*n_particle,(radius(ii)));
    
    for nz = 1:length(z)
        
        [A,B] = translate_z(max(Nmax),z(nz));
        a2 = ( A*a + B*b );
        b2 = ( A*b + B*a );
        
        pq = T * [ a2; b2 ];
        pq1 = T1 * [ a2; b2 ];
        
        p = pq(1:length(pq)/2);
        q = pq(length(pq)/2+1:end);
        
        p1 = pq1(1:length(pq)/2);
        q1 = pq1(length(pq)/2+1:end);
        
        [~,~,fz(nz),~,~,~] = forcetorque(n,m,a2,b2,p,q);
        [~,~,fz1(nz),~,~,~] = forcetorque(n,m,a2,b2,p1,q1);
        
    end
    axialrestoringforce(ii)=(min(fz));
    axialrestoringforce1(ii)=(min(fz1));
    
    zeroindex=find(fz<0,1);
    
    if length(zeroindex)~=0
        %fit to third order polynomial the local points. (only works when dz
        %sufficiently small)
        pz=polyfit(z(max([zeroindex-2,1]):min([zeroindex+2,length(z)])),fz(max([zeroindex-2,1]):min([zeroindex+2,length(z)])),2);
        root_z=roots(pz); %find roots of 3rd order poly.
        dpz=[3*pz(1),2*pz(2),1*pz(3)]; %derivative of 3rd order poly.
        
        real_z=root_z(imag(root_z)==0); % finds real roots only.
        
        rootsofsign=polyval(dpz,real_z); %roots that are stable
        zeq=real_z(rootsofsign<0); %there is at most 1 stable root. critical roots give error.
        try
            zeq=zeq(abs(zeq-z(zeroindex))==min(abs(zeq-z(zeroindex))));
        end
    else
        zeq=[];
    end
    
    if length(zeq)==0
        warning('No axial equilibrium in range!')
        zeq=0;
    end
    % equilibrium probably only correct to 1 part in 1000.
    %now work out spherical coordinates along that axis:
    [rt,theta,phi]=xyz2rtp(r,0,zeq);
    
    %calculate the x-axis coefficients for force calculation.
    Rx = rotation_matrix([-sin(0),cos(0),0],pi/2);
    Dx = wigner_rotation_matrix(max(Nmax),Rx);
    
    for nr = 1:length(r)
        
        R = rotation_matrix([-sin(phi(nr)),cos(phi(nr)),0],theta(nr)); %calculates an appropriate axis rotation off z.
        D = wigner_rotation_matrix(max(Nmax),R);
        
        [A,B] = translate_z(max(Nmax),rt(nr));
        a2 = D'*(  A * D*a +  B * D*b ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
        b2 = D'*(  A * D*b +  B * D*a ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
        
        pq = T * [ a2; b2 ];
        pq1 = T1 * [ a2; b2 ];
        p = pq(1:length(pq)/2);
        q = pq(length(pq)/2+1:end);
        p1 = pq1(1:length(pq)/2);
        q1 = pq1(length(pq)/2+1:end);
        
        [~,~,fr(nr),~,~,~] = forcetorque(n,m,Dx*a2,Dx*b2,Dx*p,Dx*q); %Dx makes the z-force calculation the x-force calculation.
        [~,~,fr1(nr),~,~,~] = forcetorque(n,m,Dx*a2,Dx*b2,Dx*p1,Dx*q1); %Dx makes the z-force calculation the x-force calculation.
        
    end
    
    transverserestoringforce(ii)=(min(fr));
    transverserestoringforce1(ii)=(min(fr1));
    
    % %test
    % figure(3)
    % plot(z,fz,'b')
    % hold on
    % plot(z,fz1,'b:')
    % plot(r,fr,'r')
    % plot(r,fr1,'r:')
    % hold off
    
end

toc

figure(1)
plot(radius*n_medium,-axialrestoringforce);
hold on
plot(radius*n_medium,-axialrestoringforce1,':');
hold off
ylabel('axial force efficiency \it Q_z')
xlabel('core radius [{\it{\lambda_{medium}}}]')
xlim(([radius(1),radius(end)]*n_medium))
ylim([0,max(ylim)])
grid on

figure(2)
plot(radius*n_medium,-transverserestoringforce);
hold on
plot(radius*n_medium,-transverserestoringforce1,':');
hold off
ylabel('transverse force efficiency \it Q_r')
xlabel('core radius [{\it{\lambda_{medium}}}]')
xlim(([radius(1),radius(end)]*n_medium))
ylim([0,max(ylim)])
grid on
