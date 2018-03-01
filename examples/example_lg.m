% Example of calculation of force in a LG03 beam trap.
%
% Approximate of Figure 2 in Nieminen et al., Optical tweezers computational
% toolbox, Journal of Optics A 9, S196-S203 (2007)
%
% How long should this take?
% The problem scales with the square of particle size. Excluding that the
% calculation goes out to main memory from cache should scale like the
% square. For the example here it took ~20 seconds on a Core2 Duo 6600 with
% 6GB of RAM.
%
% Note: the fit of the LG03 beam to the wavefunctions is rank deficient. It
% is recommended that if rank deficiency is found that the number of basis
% functions used to fit it are increased. In most cases the error is near
% round off and can be ignored. However, these errors can become quite
% large and if that is the case, a different method to map the paraxial
% scalar field to the wavefunctions should be used.
%
% PACKAGE INFO

% Specify refractive indices
n_medium = 1.33;
n_particle = 1.59;
n_relative = n_particle/n_medium;

% If you want to give all measurements in wavelengths in the surrounding
% medium, then:
wavelength = 1;
% wavelength = wavelength0 / n_medium;
% else you can give it in any units you want. Only k times lengths matters
k = 2*pi/wavelength;

% radiuz=nmax2ka([1:20]);
%
% try; clear timetakes; end
%
% for ii=1:length(radiuz)
%     tic
%
tic
radius = 2.5;
%radius = radiuz(ii);
%radius = 1;
Nmax = ka2nmax(k*radius);

diam_microns = radius * 1.064 * 2 / n_medium;
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
w0 = lg_mode_w0( [ 0 3 ], beam_angle );

% Polarisation. [ 1 0 ] is plane-polarised along the x-axis, [ 0 1 ] is
% y-polarised, and [ 1 -i ] and [ 1 i ] are circularly polarised.
polarisation = [ 1 i ];

% Location of the focal point relative to the particle. These are the
% [ x y z ] coordinates.
beam_offset = [ 0 0 0];

[n,m,a0,b0] = bsc_pointmatch_farfield(Nmax,1,[ 0 3 w0 1 polarisation 90 beam_offset ]);
[a,b,n,m] = make_beam_vector(a0,b0,n,m);

%% Insert tmatrix here %%
T = tmatrix_mie(Nmax,k,k*n_relative,radius);
%%

z = linspace(-8,8,80);
r = linspace(-4,4,80);

fz = zeros(size(z));
fr = zeros(size(r));

%root power for nomalization to a and b individually.
pwr = sqrt(sum( abs(a).^2 + abs(b).^2 ));

%normalize total momentum of wave sum to 1. Not good for SI EM field.
a=a/pwr;
b=b/pwr;

%calculate the force along z
for nz = 1:length(z)
    
    [A,B] = translate_z(Nmax,z(nz));
    a2 = ( A*a + B*b );
    b2 = ( A*b + B*a );
    
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    fz(nz) = force_z(n,m,a2,b2,p,q);
    
end

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
Rx = z_rotation_matrix(pi/2,0);
Dx = wigner_rotation_matrix(Nmax,Rx);

for nr = 1:length(r)
    
    R = z_rotation_matrix(theta(nr),phi(nr)); %calculates an appropriate axis rotation off z.
    D = wigner_rotation_matrix(Nmax,R);
    
    [A,B] = translate_z(Nmax,rt(nr));
    a2 = D'*(  A * D*a +  B * D*b ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
    b2 = D'*(  A * D*b +  B * D*a ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
    
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    fr(nr) = force_z(n,m,Dx*a2,Dx*b2,Dx*p,Dx*q); %Dx makes the z-force calculation the x-force calculation.
    
end
%     timetakes(ii)=toc;
% end
%
% plot(log([4:length(timetakes)])/log(10),log(timetakes(4:end)-timetakes(3:end-1))/log(10))
% plot([1:length(timetakes)-1],timetakes(2:end)-timetakes(1:end-1))

toc
figure; plot(z,fz);
xlabel('{\it z} (x\lambda)');
ylabel('{\it Q_z}');
aa = axis;
hold on
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');

figure; plot(r,fr);
xlabel('{\it r} (x\lambda)');
ylabel('{\it Q_r}');
aa = axis;
hold on
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
