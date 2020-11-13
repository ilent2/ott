function tests = testBscBeam
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  index = 1.33;
  omega = 2*pi*3e8 ./ 532e-9;
  bsc = ott.bsc.Bsc([1;2;3], [4;5;6]);
  beam = ott.beam.BscBeam(bsc, 'index_medium', index, 'omega', omega);
  
  testCase.verifyEqual(beam.data, bsc, 'bsc');
  testCase.verifyEqual(beam.index_medium, index, 'index');
  testCase.verifyEqual(beam.omega, omega, 'omega');
end

function testGetData(testCase)

  bsc = ott.bsc.Bsc([1;2;3;4;5;6;7;8], zeros(8, 1));
  beam = ott.beam.BscBeam(bsc);
  
  % Original beam
  trial = ott.bsc.Bsc(beam);
  testCase.verifyEqual(trial, bsc, 'original');
  
  % Grow Nmax
  testCase.verifyError(@() ott.bsc.Bsc(beam, 10), ...
    'ott:beam:BscBeam:recalculate_not_implemented');
  
  % Rotate
  rbsc = bsc.rotateY(pi/2);
  rbeam = beam.rotateY(pi/2);
  trial = ott.bsc.Bsc(rbeam);
  testCase.verifyEqual(trial, rbsc, 'rotate');
  
  % Translate small
  tbsc = bsc.translateXyz([0.001;0;0], 'Nmax', 1);
  tbeam = beam.translateXyz([0.001;0;0]*beam.wavelength);
  trial = ott.bsc.Bsc(tbeam);
  testCase.verifyEqual(trial, tbsc, 'translate');

end

function testPlaneWaveFields(testCase)

  beam = ott.beam.PlaneWave('Nmax', 20);
  
  % Check E at two points
  E0 = beam.efield([0;0;0]).vxyz;
  E1 = beam.efield([0;0;1]*beam.wavelength).vxyz;
  testCase.verifyEqual(E0, E1, 'AbsTol', 1e-6, 'RelTol', 1e-6, ...
    'fields should be the same when separated by wavelength');
  
  % Check derivative of E
  dz = 1.0e-3*beam.wavelength;
  E1 = beam.efield([0;0;dz]).vxyz;
  E2 = beam.efield([0;0;-dz]).vxyz;
  Edd = (E1 - 2*E0 + E2)./dz^2;
  testCase.verifyEqual(abs(Edd), abs(beam.omega^2/beam.speed^2.*E0), ...
    'AbsTol', 2e-5, 'RelTol', 1e-5, 'wave equation for E not satisfied');
  
  % Check H relative to E
  H0 = beam.hfield([0;0;0]).vxyz;
  testCase.verifyEqual(abs(H0([2,1,3])), abs(E0./beam.impedance), ...
    'AbsTol', 1e-6, 'H should be scaled by beam impedance');
  
  % Check H relative to E (in the far-field)
  S = warning('off', 'ott:beam:BscInfinite:farfield_is_finite');
  testCase.addTeardown(@() warning(S));
  E0 = beam.efarfield([0;0;0]).vxyz;
  H0 = beam.hfarfield([0;0;0]).vxyz;
  testCase.verifyEqual(abs(H0([2,1,3])), abs(E0./beam.impedance), ...
    'AbsTol', 1e-6, 'FF: H should be scaled by beam impedance');

end

function testNearfieldFunctionsCoverage(testCase)

  beam = ott.beam.PlaneWave('Nmax', 20);
  
  xyz = randn(3, 5)*beam.wavelength;
  rtp = ott.utils.xyz2rtp(xyz);
  
  trialX = beam.efieldRtp(rtp);
  trialR = beam.efield(xyz);
  testCase.verifyEqual(trialX.vrtp, trialR.vrtp, 'RelTol', 1e-15, 'efield');
  
  trialX = beam.hfieldRtp(rtp);
  trialR = beam.hfield(xyz);
  testCase.verifyEqual(trialX.vrtp, trialR.vrtp, 'RelTol', 1e-15, 'hfield');
  
  [trialEX, trialHX] = beam.ehfieldRtp(rtp);
  [trialER, trialHR] = beam.ehfield(xyz);
  testCase.verifyEqual(trialEX.vrtp, trialER.vrtp, 'RelTol', 1e-15, 'ehfield e');
  testCase.verifyEqual(trialHX.vrtp, trialHR.vrtp, 'RelTol', 1e-15, 'ehfield h');
  
end

function testFarfieldFunctions(testCase)

  index = 1.33;
  omega = 2*pi*3e8 ./ 532e-9;
  bsc = ott.bsc.Bsc([1;2;3], [4;5;6]);
  beam = ott.beam.BscBeam(bsc, 'index_medium', index, 'omega', omega);
  
  rtp = randn(3, 5);
  
  trialE = beam.efarfield(rtp);
  trialH = beam.hfarfield(rtp);
  
  [trialEX, trialHX] = beam.ehfarfield(rtp);
  testCase.verifyEqual(trialEX.vrtp, trialE.vrtp, 'RelTol', 1e-15, 'eh e');
  testCase.verifyEqual(trialHX.vrtp, trialH.vrtp, 'RelTol', 1e-15, 'eh h');
  
end

function testForceFunctions(testCase)
  % These are effectively definitions from OTTv1

  bsc1 = ott.bsc.Bsc([1;2;3], 0);
  bsc2 = ott.bsc.Bsc([3;2;1], 0);
  
  beam1 = ott.beam.BscBeam(bsc1);
  beam2 = ott.beam.BscBeam(bsc2);
  
  target = bsc1.force(bsc2) ./ beam1.speed;
  trial = beam1.force(beam2);
  testCase.verifyEqual(target, trial, 'force');
  
  target = bsc1.torque(bsc2) ./ beam1.omega;
  trial = beam1.torque(beam2);
  testCase.verifyEqual(target, trial, 'torque');
  
  target = bsc1.spin(bsc2) ./ beam1.omega;
  trial = beam1.spin(beam2);
  testCase.verifyEqual(target, trial, 'spin');

end

function testScatter(testCase)

  n_medium = 1.35;
  n_rel = 1.2;
  beam = ott.beam.PlaneWave('index_medium', n_medium);
  
  sphere = ott.shape.Sphere(beam.wavelength);
  particle = ott.particle.Fixed.FromShape(sphere, ...
    'index_relative', n_rel, 'internal', true, ...
    'wavelength0', beam.wavelength0);
  
  sbeam = beam.scatter(particle);
  testCase.assertInstanceOf(sbeam, 'ott.beam.Scattered');
  testCase.verifyEqual(sbeam.index_medium, n_medium, 'sbeam_index');
  testCase.verifyEqual(sbeam.scattered.index_medium, n_medium, 'scat index');
  testCase.verifyEqual(sbeam.incident.index_medium, n_medium, 'inc index');
  testCase.verifyEqual(sbeam.internal.index_medium, n_medium * n_rel, 'int index');

end

function testScatterForce(testCase)

  beam = ott.beam.Gaussian();
  particle = ott.particle.Fixed.FromShape(...
    ott.shape.Sphere(beam.wavelength0), 'index_relative', 1.2, ...
    'wavelength0', beam.wavelength0);
  
  pos = [0.5;0;0]*beam.wavelength;
  
  sbeam = beam.scatter(particle, 'position', pos);
  F1 = sbeam.force();
  
  F2 = beam.force(particle, 'position', pos);
  
  % Need to swap position (get force in particle reference fraome)
  beam.position = -sbeam.position;
  sbeam.position = [0;0;0];
  F3 = (sbeam.intensityMoment('type', 'scattered') ...
    - beam.intensityMoment())./beam.speed;
  
  testCase.verifyEqual(F1, F2, 'RelTol', 6e-2, ...
    'direct force and scatter force dont match');
  testCase.verifyEqual(F3, F2, 'RelTol', 6e-2, ...
    'farfield integral doesn''t match');

end

function testPowerIntegralGaussianNearfield(testCase)

  waist = 1e-6;
  beam = ott.beam.Gaussian(waist, ...
    'truncation_angle', 0.8*pi/2);
  
  W = beam.wavelength*200;
  power = integral(@calcS, 0, W, 'RelTol', 1e-3);
  
  testCase.verifyEqual(power, beam.power, 'AbsTol', 1e-2);

  function s = calcS(x)
    [E, H] = beam.ehfield([x(:),0*x(:),0*x(:)].');
    S = real(cross(E, conj(H)));
    s = reshape(sum([0;0;1].*S, 1), size(x));
    s = s .* 2*pi.*x;
  end
end


