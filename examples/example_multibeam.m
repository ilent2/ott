% Calculates the force from multiple (identical) beams positioned at
% different locations.  Genrates a figure which compares forces calculated
% using incoherent beams and two methods of adding coherent beams.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Make warnings less obtrusive
ott_warning('once');
change_warnings('off');

n_relative = 1.2; %relative refractive index = 1.2
r_particle = 1;   %particle radius is 1 wavelength in the medium
beam_angle = 50;  %beam half angle of 50 degrees
polarization = [1 i];  %circular polarization

Nmax = ka2nmax(2*pi*r_particle); %calculate limit of the expansion

%calculate the beam shape coefficients
[n, m, a0, b0] = bsc_pointmatch_farfield(Nmax, 1, [0 0 beam_angle 1 polarization 90 ]);
[a, b, n, m] = make_beam_vector(a0, b0, n, m); %pack beam vector to full size

%calculate the power of the two beam system:
power_total = 2 * sum(abs(a).^2 + abs(b).^2);

%% Coherent beams: method 1 - Re-size Nmax, separate beams 1lambda in the x-direction
%Let's make the a's and b's bigger by 1 wavelength (kr~12) each.
Nmax=Nmax+12;
aa = change_nmax(a,Nmax);
bb = change_nmax(b,Nmax);

%these beams have the same properties, but are phase shifted by pi:
a1 = aa;
b1 = bb;
a2 = aa * exp(2 * pi * 1i * 1/2); %shift half a wavelength (in z) in the medium.
b2 = bb * exp(2 * pi * 1i * 1/2); %these phase shifts are for demonstration only...

%calculate a rotation off the +z-axis
R = rotation_matrix([-sin(0),cos(0),0],pi/2);
D = wigner_rotation_matrix(Nmax,R);

%shift the beams so that they are symmetric around the origin:
[A,B] = translate_z(Nmax,1);

%rotate, translate and rotate' beam 1
a10 = D' * (A * D * a1 + B * D * b1);
b10 = D' * (A * D * b1 + B * D * a1);

%move beam 2. Either reverse the rotation or translation.
a20 = D * (A * D' * a2 + B * D' * b2);
b20 = D * (A * D' * b2 + B * D' * a2);

%calculate the total field BSC:
a_total = a10 + a20;
b_total = b10 + b20;

T = tmatrix_mie(Nmax,2*pi,2*pi*n_relative,r_particle);   %T-matrix for a Mie scatterer

pq = T * [a_total; b_total];                %calculate scattered BSC

%our n's and m's are the wrong size... Let's recalculate them!
%There are nmax^2 + 2 * nmax modes for a given truncation, so:
[n1,m1]=combined_index([1:Nmax^2 + 2 * Nmax].');

%calculate force at point in z.
[~,~,Q_z,~,~,~] = forcetorque(n1, m1, [a_total; b_total], pq);
Q_z = Q_z / power_total;

%demonstrate success by plotting a line of force along the separation direction.
figure(1)
clf
hold on
for ii=1:100
    [A,B] = translate_z(Nmax,4*ii/100-2); %move from [-2,2]
    % rotate, translate and rotate' the beam so we can calculate the
    % scattering as a function of position
    a_n = D' * (A * D * a_total + B * D * b_total); 
    b_n = D' * (A * D * b_total + B * D * a_total);
    
    pq=T*[a_n;b_n];
    
    %plot the force in the x direction by rotating the spherical waves
    [~,~,Q_x,~,~,~] = forcetorque(n1, m1, D*a_n, D*b_n, ...
        D*pq(1:end/2), D*pq(end/2+1:end));
    Q_x = Q_x / power_total;
    plot(4*ii/100-2,Q_x,'bx')
end
hold off
vech=get(gca,'children'); %get the handles on the axis.
h(1)=vech(1);             %extract a handle for the legend.

%% Coherent beams: method 2 - Separate beams 1lambda in the x-direction using only local regions.
Nmax = ka2nmax(2*pi*r_particle); %re-calculate limit of the expansion
T = tmatrix_mie(Nmax,2*pi,2*pi*n_relative,r_particle);   %T-matrix for a Mie scatterer

%calculate a rotation off the +z-axis
R = rotation_matrix([-sin(0),cos(0),0],pi/2);
D = wigner_rotation_matrix(Nmax,R);

av=[a,a]; %using old nmax
bv=[b,b]; %ditto
beam_locations=[-1,1;0,0;0,0]; %beam location in [x;y;z].

%unlike the previous example the important part where we add the beams is
%still to come
figure(1)
hold on
for ii=1:100
    %this method requires a second loop inside the first, this is where the
    %total scattering of the sphere is calculated.
    for jj=1:size(beam_locations,2)
        % Calculate the translation to get the particle to the correct
        % location in space relative to the beam.
        [A,B] = translate_z(Nmax,4*ii/100-2-beam_locations(1,jj));
        if jj==1 %the first loop starts with a beam with no power
            a_n = D' * (A * D * av(:,1) + B * D * bv(:,1));
            b_n = D' * (A * D * bv(:,1) + B * D * av(:,1));
        else     %each succesive loop adds a beam
            a_n = a_n + D' * (A * D * av(:,jj) + B * D * bv(:,jj));
            b_n = b_n + D' * (A * D * bv(:,jj) + B * D * av(:,jj));
        end
    end
    
    %calculate the scattering...
    pq=T*[a_n;b_n];
    
    [~,~,Q_x,~,~,~] = forcetorque(n, m, D*a_n, D*b_n, ...
        D*pq(1:end/2), D*pq(end/2+1:end));
    Q_x = Q_x / power_total;
    plot(4*ii/100-2,Q_x,'rx')
end
hold off
vech=get(gca,'children'); %get the handles on the axis.
h(2)=vech(1);             %extract a handle for the legend.

%% Incoherent beams: adding forces...
Nmax = ka2nmax(2*pi*r_particle); %re-calculate limit of the expansion
T = tmatrix_mie(Nmax,2*pi,2*pi*n_relative,r_particle);   %T-matrix for a Mie scatterer

%calculate a rotation off the +z-axis
R = rotation_matrix([-sin(0),cos(0),0],pi/2);
D = wigner_rotation_matrix(Nmax,R);

av=[a,a]; %using old nmax
bv=[b,b]; %ditto
beam_locations=[-1,1;0,0;0,0]; %beam location in [x;y;z].

figure(1)
hold on
for ii=1:100
    %this method requires a second loop inside the first, this is where the
    %total force of the sphere is calculated.
    for jj=1:size(beam_locations,2)
        [A,B] = translate_z(Nmax,4*ii/100-2-beam_locations(1,jj));
        
        a_n = D' * (A * D * av(:,jj) + B * D * bv(:,jj));
        b_n = D' * (A * D * bv(:,jj) + B * D * av(:,jj));
        
        pq=T*[a_n;b_n];
        
        if jj==1
            [~,~,Q_x,~,~,~] = forcetorque(n, m, D*a_n, D*b_n, ...
                D*pq(1:end/2), D*pq(end/2+1:end));
            Q_x = Q_x / power_total;
        else
            [~,~,temp_x,~,~,~]=forcetorque(n, m, D*a_n, D*b_n, ...
                D*pq(1:end/2), D*pq(end/2+1:end));
            Q_x = Q_x + temp_x / power_total;
        end
    end
    
    plot(4*ii/100-2,Q_x,'ro')
end
hold off
vech=get(gca,'children');
h(3)=vech(1);

%% Incoherent beams: demonstrate that phase has no effect on incoherent beam method
av=[a,a*exp(2*pi*1i*1/2)]; %phase shift the second beam wrt the first
bv=[b,b*exp(2*pi*1i*1/2)]; %phase shift the second beam wrt the first

figure(1)
hold on
for ii=1:100
    %this method requires a second loop inside the first, this is where the
    %total force of the sphere is calculated.
    for jj=1:size(beam_locations,2)
        [A,B] = translate_z(Nmax,4*ii/100-2-beam_locations(1,jj));
        
        a_n = D' * (A * D * av(:,jj) + B * D * bv(:,jj));
        b_n = D' * (A * D * bv(:,jj) + B * D * av(:,jj));
        
        pq=T*[a_n;b_n];
        
        if jj==1
            [~,~,Q_x,~,~,~] = forcetorque(n, m, D*a_n, D*b_n, ...
                D*pq(1:end/2), D*pq(end/2+1:end));
            Q_x = Q_x / power_total;
        else
            [~,~,temp_x,~,~,~]=forcetorque(n, m, D*a_n, D*b_n, ...
                D*pq(1:end/2), D*pq(end/2+1:end));
            Q_x = Q_x + temp_x / power_total;
        end
    end
    
    plot(4*ii/100-2,Q_x,'bo')
end
hold off
vech=get(gca,'children');
h(4)=vech(1);

%% Add legends and axis labels to figure

grid on
xlabel('x [\lambda]')
ylabel('Q_x')
M={'large N_{max} coherent (phase shift)', ...
    'small N_{max} coherent (no shift)', ...
    'incoherent (no shift)','incoherent (phase shift)'};
legend(h,M);
axis tight
