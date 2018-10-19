function tests = examples
  tests = functiontests(localfunctions);
end

function testSphere(testCase)

  addpath('../');
  ott.change_warnings('off');

  % Refractive index of particle and medium
  n_medium = 1.33;
  n_particle = 1.59;

  % Wavelength of light in vacuum [m]
  wavelength0 = 1064e-9;

  % Calculate the wavelength in the medium
  wavelength_medium = wavelength0/n_medium;

  % Radius of particle
  radius = 1.0*wavelength_medium;

  % Specify the numerical aparture of the beam/microscope
  NA = 1.02;

  % Create a T-matrix for a sphere
  T = ott.Tmatrix.simple('sphere', radius, 'wavelength0', wavelength0, ...
      'index_medium', n_medium, 'index_particle', n_particle);

  % Create a simple Gaussian beam with circular polarisation
  beams = {};
  beams{1} = ott.BscPmGauss('NA', NA, 'polarisation', [ 1 1i ], ...
      'index_medium', n_medium, 'wavelength0', wavelength0);

  beam_angle = asin(NA/n_medium);
  beams{2} = ott.BscPmGauss('lg', [ 0 3 ], ...
      'polarisation', [ 1 1i ], 'angle', beam_angle, ...
      'index_medium', n_medium, 'wavelength0', wavelength0);

  beam_angle = asin(NA/n_medium);
  beams{3} = ott.BscPmGauss('hg', [ 2 3 ], ...
      'polarisation', [ 1 1i ], 'angle', beam_angle, ...
      'index_medium', n_medium, 'wavelength0', wavelength0);

  npts = 20;
  warning_status = { false, false, true };
  tol = 0.01;
  fr_target = {};
  fr_target{1} = [0.000044433787344 0.000083238173655 0.000273170525965, 0.000730342126549, ...
    0.001386727616596   0.006564585296498   0.028535350073358,  0.131448283401034, ...
    0.067949468062682   0.001883673983584   0.001883673983584, 0.067949468062682, ...
    0.131448283401034   0.028535350073358   0.006564585296498, 0.001386727616596, ...
    0.000730342126549   0.000273170525965   0.000083238173655   0.000044433787344 ];
  fr_target{2} = [0.004701026373629   0.008546824420202   0.014118834900474   0.020515811480778, ...
    0.023190400451052   0.015981070354735   0.000737396239011  -0.009038328405627, ...
    -0.010507385325924  -0.003761950950816  -0.003761950950816  -0.010507385325924, ...
    -0.009038328405627   0.000737396239011   0.015981070354735   0.023190400451052, ...
    0.020515811480778   0.014118834900474   0.008546824420203   0.004701026373629 ];
  fr_target{3} = [0.000010117098003   0.000041306173553   0.000091951242659   0.000910668824001, ...
    0.005802261157539   0.007680208357796   0.017103312085298 0.028525744768119, ...
    0.019069772279068   0.009485957013539   0.009485957013539   0.019069772279068, ...
    0.028525744768119   0.017103312085298   0.007680208357793   0.005802261157539, ...
    0.000910668824001   0.000091951242659   0.000041306173553   0.000010117098003 ];

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.RelativeTolerance;

  for ii = 1:length(beams)

    beam = beams{ii};

    % Normalise power
    beam.power = 1.0;

    %calculate the force along z
    z = [0;0;1]*linspace(-8,8,npts)*wavelength_medium;
    fz = ott.forcetorque(beam, T, 'position', z);

    % Find the equilibrium along the z axis
    zeq = ott.find_equilibrium(z(3, :), fz(3, :));
    testCase.verifyThat(isempty(zeq), IsEqualTo(warning_status{ii}));
    if isempty(zeq)
      zeq=0;
    end
    zeq = zeq(1);

    % Calculate force along x-axis (with z = zeq, if found)
    r = [1;0;0]*linspace(-4,4,npts)*wavelength_medium + [0;0;zeq];
    fr = ott.forcetorque(beam, T, 'position', r);

    testCase.verifyThat(fr(3, :), IsEqualTo(fr_target{ii}, ...
        'Within', RelativeTolerance(tol)));
  end

end

function testAxialEquilibrium(testCase)

  addpath('../');
  ott.change_warnings('off');

  % Specify refractive indices
  n_medium = 1.34;
  n_particle = 1.59;

  % Specify the wavelength in freespace [m]
  wavelength = 1064.0e-9;

  % Specify the particle radius (sphere)
  radius = 1.5*wavelength/n_medium;

  %% Calculate the beam

  beam = ott.BscPmGauss('angle_deg', 50, 'polarisation', [ 1 0 ], ...
      'index_medium', n_medium, 'wavelength0', wavelength, 'power', 1.0);

  %% Calculate the particle T-matrix

  T = ott.Tmatrix.simple('sphere', radius, ...
      'index_medium', n_medium, ...
      'index_particle', n_particle, ...
      'wavelength0', wavelength);

  %% Find the equilibrium and trap stiffness in the x and x directions

  % Find the equilibrium in the z-direction
  [z,kz] = ott.axial_equilibrium(T, beam);

  % Translate the beam to the z-axis equilibrium
  beam = beam.translateZ(z);

  % Rotate the beam about the y axis (so the beam is aligned with the x axis)
  beam = beam.rotateY(pi/2.0);

  % Calculate the equilibrium in the x-direction
  [x,kx] = ott.axial_equilibrium(T, beam);

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.RelativeTolerance;
  import matlab.unittest.constraints.AbsoluteTolerance;
  
  tol = 0.01;
  tolm = 1e-9;

  testCase.verifyThat(x, IsEqualTo(0.0, ...
        'Within', AbsoluteTolerance(tolm)));

  testCase.verifyThat(kx, IsEqualTo(-2.450214062064871e+05, ...
        'Within', RelativeTolerance(tol)));

  testCase.verifyThat(z, IsEqualTo(4.023474767793779e-07, ...
        'Within', AbsoluteTolerance(tolm)));

  testCase.verifyThat(kz, IsEqualTo(-6.258368435930368e+04, ...
        'Within', RelativeTolerance(tol)));
end

function testCube(testCase)

  % Specify refractive indices
  n_medium = 1.33;
  n_particle = 1.57;

  % Specify the wavelength in the freespace [m]
  wavelength0 = 1064e-9;
  wavelength_medium = wavelength0 / n_medium;

  % Specify the size of the cube
  width = 1.0*wavelength_medium;

  beam = ott.BscPmGauss('NA', 1.25, 'polarisation', [ 1i 1 ], ...
    'power', 1.0, 'index_medium', n_medium, 'wavelength0', wavelength0);

  T = ott.Tmatrix.simple('cube', width, ...
        'index_medium', n_medium, ...
        'index_particle', n_particle, ...
        'wavelength0', wavelength0);

  x = [ 0; 1; 1 ];
  Rw = rotz(45)*roty(45);

  %% Calculate the force and torque by matching Nmax

  % Ensure Nmax of beam and particle are the same (optimisation)
  T1 = T;
  beam1 = beam;
  Nmax = max(beam.Nmax, T1.Nmax);
  T1.Nmax = Nmax;
  beam1.Nmax = Nmax;

  % Translate the beam to the particle
  tbeam = beam1.translateXyz(x * wavelength_medium);

  % Rotate the beam to match the particle rotation
  [rbeam, D] = tbeam.rotate(Rw);

  % Scatter the beam and return from the particle frame to the lab frame
  sbeam = D' * (T1 * rbeam);

  % Calculate the force and the torque on the particle
  [ff1,tt1] = ott.forcetorque(tbeam, sbeam);

  %% Calculate the force and torque by shrinking the beam Nmax

  % Translate the beam to the particle
  tbeam = beam.translateXyz(x * wavelength_medium, 'Nmax', T.Nmax);

  % Rotate the beam to match the particle rotation
  [rbeam, D] = tbeam.rotate(Rw);

  % Scatter the beam and return from the particle frame to the lab frame
  sbeam = D' * (T * rbeam);

  % Calculate the force and the torque on the particle
  [ff2,tt2] = ott.forcetorque(tbeam, sbeam);

  %% Calculate the force and torque in one line

  [ff3,tt3] = ott.forcetorque(beam, T, ...
      'position', x * wavelength_medium, 'rotation', Rw);

  %% Check the calculated forces/torques match

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.RelativeTolerance;
  tol = 1.0e-3;

  testCase.verifyThat(ff1, IsEqualTo(ff2, ...
        'Within', RelativeTolerance(tol)), 'Multiline force');

  testCase.verifyThat(tt1, IsEqualTo(tt2, ...
        'Within', RelativeTolerance(tol)), 'Multiline torque');

  testCase.verifyThat(tt3, IsEqualTo(tt2, ...
        'Within', RelativeTolerance(tol)), 'One line force');

  testCase.verifyThat(tt3, IsEqualTo(tt2, ...
        'Within', RelativeTolerance(tol)), 'One line torque');

end
