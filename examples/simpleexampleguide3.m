%code to calculate counter propagating beams
n_relative = 2; %relative refractive index = 2
r_particle = 1;   %particle radius is 1 wavelength in the medium
beam_angle = 50;  %beam half angle of 50 degrees
polarization = [1 0];  %circular polarization
beam_offset = [0 0 0]; %sets the destination of the beam focus to O.

Nmax = ka2nmax(2*pi*r_particle); %calculate limit of the expansion

w0 = lg_mode_w0([0 0], beam_angle);

%calculate the beam shape coefficients
[n, m, a0, b0] = bsc_pointmatch_farfield(Nmax, 1, [0 0 w0 1 polarization 90 beam_offset]);
[a, b, n, m] = make_beam_vector(a0, b0, n, m); %pack beam vector to full size

%calculate the power of the two beam system:
power_total = 2 * sum(abs(a).^2 + abs(b).^2);

%we want to generate a system of counter propagaing beams. To do this we will need to rotate the beam.
R=z_rotation_matrix(pi, 0);
D=wigner_rotation_matrix(Nmax, R);

%calculate the two beams.
a1=a;
b1=b;
a2=D*a;
b2=D*b;

%note that the rotation we did also flips x because of the way we have defined the rotation.
a_total=a1+a2;
b_total=b1+b2;

T = tmatrix_mie(Nmax,2*pi,2*pi*n_relative,r_particle);   %T-matrix for a Mie scatterer

Q_z = zeros(100,1);
z = linspace(-2,2,100).';
figure(1)
clf;
for ii=1:100
    [A,B]=translate_z(Nmax,z(ii));
    a_n=A*a_total+B*b_total;
    b_n=A*b_total+B*a_total;
    
    pq = T * [a_n; b_n];
    
    Q_z(ii) = force_z(n, m, [a_n; b_n], pq) / power_total;
end
plot(z,Q_z,'b')

%make a vector of vectors for a and b
av=[a1,a2];
bv=[b1,b2];

Q_z1 = zeros(100,2);
z = linspace(-2,2,100).';
for ii=1:100
    for jj=1:2
        [A,B]=translate_z(Nmax,z(ii));
        a_n=A*av(:,jj)+B*bv(:,jj);
        b_n=A*bv(:,jj)+B*av(:,jj);
        
        pq = T * [a_n; b_n];
        
        Q_z1(ii,jj) = force_z(n, m, [a_n; b_n], pq) / power_total;
    end
end

figure(1)
hold on;
plot(z,sum(Q_z1,2),'rx')
plot(z,(Q_z1(:,1)),'go')
plot(z,(Q_z1(:,2)),'mo')
grid on
hold off
xlabel('z [\lambda]','fontsize',16)
ylabel('Q_z','fontsize',16)
set(gca,'fontsize',16)
legend('counter propagating','incoherent counter propagating')