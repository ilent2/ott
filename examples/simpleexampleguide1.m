%code to calculate the force acting upon a wavelength radius sphere at the origin
n_relative=1.2; %relative refractive index = 1.2
r_particle=1;   %particle radius is 1 wavelength in the medium
beam_angle=50;  %beam half angle of 50 degrees
polarization = [1 i];  %circular polarization
beam_offset = [0 0 0]; %sets the destination of the beam focus to O.

Nmax = ka2nmax(2*pi*r_particle); %calculate limit of the expansion

w0 = lg_mode_w0([0 0], beam_angle);

%calculate the beam shape coefficients
[n, m, a0, b0] = bsc_pointmatch_farfield(Nmax, 1, [0 0 w0 1 polarization 90 beam_offset]);
[a, b, n, m] = make_beam_vector(a0, b0, n, m); %pack beam vector to full size

T = tmatrix_mie(Nmax,2*pi,2*pi*n_relative,r_particle);   %T-matrix for a Mie scatterer

pq = T * [a; b];                %calculate scattered BSC

%finally, calculate the force in the xyz directions
Q_xyz = forcetorque(n, m, [a; b], pq) / sum(abs(a).^2 + abs(b).^2);

%now calculate the transformations which move the particle one wavelength
[A, B] = translate_z(Nmax,1);

a_new = A * a + B * b;
b_new = A * b + B * a;

pq_new = T * [a_new; b_new];

%calculate the force on the z-axis
Q_z = force_z(n, m, [a_new; b_new], pq_new) / sum(abs(a).^2 + abs(b).^2);

%calculate a rotation off the +z-axis
R = z_rotation_matrix(pi/4,pi/4);
D = wigner_rotation_matrix(Nmax,R);

%rotate the BSC
a_r = D * a;
b_r = D * b;

%translate the BSC in the defined direction
a_rtemp = A * a_r + B * b_r;
b_rtemp = A * b_r + B * a_r;

%rotate back to the original coordinate rotation.
a_rnew = D'*a_rtemp;
b_rnew = D'*b_rtemp;

pq_rnew = T * [a_rnew; b_rnew];

Q_xyz = forcetorque(n, m, [a_rnew; b_rnew], pq_rnew) / sum(abs(a).^2 + abs(b).^2);
