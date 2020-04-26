function tests = forcetorque
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)

  % Ensure the ott package is in our path
  addpath('../');

  % Generate a gaussian beam for testing
  testCase.TestData.beam = ott.BscPmGauss('power', 1.0);

  % Generate a spherical particle for testing
  testCase.TestData.T = ott.Tmatrix.simple('sphere', 1.0, ...
      'index_relative', 1.2);

  % Non-contrasting particle
  testCase.TestData.Tnocont = ott.Tmatrix.simple('sphere', 1.0, ...
      'index_relative', 1.0);

end

function testSimple(testCase)

  T = testCase.TestData.T;
  beam = testCase.TestData.beam;

  % Calculate scattered beam
  sbeam = T * beam;

  % Calculate force and torque
  [f, t] = ott.forcetorque(beam, sbeam);

  % TODO: Check it works

end

function testUnevenNmax(testCase)

  beam = testCase.TestData.beam;
  T = testCase.TestData.T;
  T.Nmax = [beam.Nmax+2, beam.Nmax+1];

  % Calculate force and torque
  [f0, t0] = ott.forcetorque(beam, testCase.TestData.T);
  [f1, t1] = ott.forcetorque(beam, T);
  [f2, t2] = ott.forcetorque(beam, [T, T]);

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  testCase.verifyThat(f1, IsEqualTo(f0, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T1 force does not match');

  testCase.verifyThat(f2(1:3), IsEqualTo(f0, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T2(1) force does not match');
  testCase.verifyThat(f2(4:6), IsEqualTo(f0, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T2(2) force does not match');

end

function testLargeParticles(testCase)

  beam = testCase.TestData.beam;
  T = testCase.TestData.T;
  T.Nmax = beam.Nmax+2;

  % Calculate force and torque
  [f0, t0] = ott.forcetorque(beam, testCase.TestData.T);
  [f1, t1] = ott.forcetorque(beam, T);
  [f2, t2] = ott.forcetorque(beam, [T, T]);

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  testCase.verifyThat(f1, IsEqualTo(f0, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T1 force does not match');

  testCase.verifyThat(f2(1:3), IsEqualTo(f0, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T2(1) force does not match');
  testCase.verifyThat(f2(4:6), IsEqualTo(f0, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T2(2) force does not match');
end

function testMultiParticle(testCase)

  beam = testCase.TestData.beam;
  T = testCase.TestData.T;
  T1 = ott.Tmatrix.simple('sphere', 0.8, ...
      'index_relative', 1.2);

  % Calculate force and torque
  [f, t] = ott.forcetorque(beam, [T, T1]);

  % Compare the result to the individual beams (requires other tests to pass)
  [f1, t1] = ott.forcetorque(beam, T);
  [f2, t2] = ott.forcetorque(beam, T1);

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  testCase.verifyThat(f1, IsEqualTo(f(1:3), ...
      'Within', AbsoluteTolerance(tol)), ...
      'T1 force does not match');

  testCase.verifyThat(f2, IsEqualTo(f(4:6), ...
      'Within', AbsoluteTolerance(tol)), ...
      'T2 force does not match');

end

function testPlaneWave(testCase)

  Nmax = 12;

  % Test the force on a non-contrasting particle
  T = testCase.TestData.Tnocont;
  beam = ott.BscPlane(0.0, 0.0, 'Nmax', Nmax, 'polarisation', [1, -1i]);
  beam.power = 1.0;
  
  numpts = 10;
  R = zeros(3, 3*numpts);
  for ii = 1:numpts
    R(:, (1:3) + 3*(ii-1)) = rotz(rand()*360)*roty(rand()*180);
  end
  
  [f, ~] = ott.forcetorque(beam, T, 'rotation', R);
  
  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-3;

  testCase.verifyThat(vecnorm(f), IsEqualTo(zeros(1, numpts), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Force does not vanish for non-contrasing particle');
    
  % Create an identity T-matrix in the total field type
  T = ott.Tmatrix(-0.5*speye(2*ott.utils.combined_index(Nmax, Nmax)), ...
      'scattered');
  
  [f, t, s] = ott.forcetorque(beam, T, 'rotation', R);

  testCase.verifyThat(vecnorm(f), IsEqualTo(0.9231*ones(1, numpts), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Force should approach 0.9231 in this configuration');

  testCase.verifyThat(vecnorm(t), IsEqualTo(ones(1, numpts), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Torque should approach 1 in this configuration');

  testCase.verifyThat(vecnorm(s), IsEqualTo(0.9231*ones(1, numpts), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Spin should approach 0.9231 in this configuration');

end

function testMultiBeam(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  T = testCase.TestData.T;
  beam = testCase.TestData.beam;

  % Add a second beam
  pbeam = beam.append(beam.translateZ(pi/2));

  testCase.verifyThat(pbeam.Nbeams, IsEqualTo(2), ...
      'Nbeams after packing failed');
  testCase.verifyThat(pbeam.Nmax, IsEqualTo(beam.Nmax), ...
      'Nmax after packing failed');

  % Calculate scattered beam
  sbeam = T * pbeam;

  testCase.verifyThat(sbeam.Nbeams, IsEqualTo(2), ...
      'Beams after scattering incorrect');
  testCase.verifyThat(sbeam.Nmax, IsEqualTo(min(T.Nmax(1), beam.Nmax)), ...
      'Nmax after scattering incorrect');

  % Calculate force and torque
  [f, t] = ott.forcetorque(beam, T, 'position', [0;0;1]*[0, pi/2]);

  % Compare the result to the individual beams (requires other tests to pass)
  [f1, t1] = ott.forcetorque(pbeam.beam(1), sbeam.beam(1));
  [f2, t2] = ott.forcetorque(pbeam.beam(2), sbeam.beam(2));

  testCase.verifyThat(size(f), IsEqualTo([3, 2]), ...
      'Beam1 force size incorrect (1)');

  testCase.verifyThat(f1(:), IsEqualTo(f(1:3).', ...
      'Within', AbsoluteTolerance(tol)), ...
      'Beam1 force does not match (1)');

  testCase.verifyThat(f2(:), IsEqualTo(f(4:6).', ...
      'Within', AbsoluteTolerance(tol)), ...
      'Beam2 force does not match (1)');

  % Calculate force and torque from two arrays of beams
  [f, t] = ott.forcetorque(pbeam, sbeam);

  testCase.verifyThat(size(f), IsEqualTo([3, 2]), ...
      'Beam1 force size incorrect (2)');

  testCase.verifyThat(f1(:), IsEqualTo(f(1:3).', ...
      'Within', AbsoluteTolerance(tol)), ...
      'Beam1 force does not match (2)');

  testCase.verifyThat(f2(:), IsEqualTo(f(4:6).', ...
      'Within', AbsoluteTolerance(tol)), ...
      'Beam2 force does not match (2)');

  % Calculate force and torque from T-matrix and array of beams
  [f, t] = ott.forcetorque(pbeam, T, 'position', [0;0;1]*[0, pi/2]);

  beam3 = pbeam.beam(1).translateZ(pi/2);
  beam4 = pbeam.beam(2).translateZ(pi/2);
  [f3, t1] = ott.forcetorque(beam3, T*beam3);
  [f4, t2] = ott.forcetorque(beam4, T*beam4);

  testCase.verifyThat(size(f), IsEqualTo([3, 2, 2]), ...
      'Beam1 force size incorrect (2)');

  testCase.verifyThat(f1(:), IsEqualTo(f(1:3).', ...
      'Within', AbsoluteTolerance(tol)), ...
      'Beam1 force does not match (3)');

  testCase.verifyThat(f2(:), IsEqualTo(f(4:6).', ...
      'Within', AbsoluteTolerance(tol)), ...
      'Beam2 force does not match (3)');

  testCase.verifyThat(f3(:), IsEqualTo(f(7:9).', ...
      'Within', AbsoluteTolerance(tol)), ...
      'Beam3 force does not match (3)');

  testCase.verifyThat(f4(:), IsEqualTo(f(10:12).', ...
      'Within', AbsoluteTolerance(tol)), ...
      'Beam4 force does not match (3)');
end

function testCoherent(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  T = testCase.TestData.T;
  beam = testCase.TestData.beam;

  % Add a second beam
  pbeam = beam.append(beam.translateZ(pi/2));

  [f, t] = ott.forcetorque(pbeam, T, ...
    'position', [0;0;1]*[0, pi/2], 'coherent', true);
  
  beam1 = pbeam.mergeBeams();
  beam2 = pbeam.translateZ(pi/2).mergeBeams();
  [f1, t1] = ott.forcetorque(beam1, T);
  [f2, t2] = ott.forcetorque(beam2, T);
  
  testCase.verifyThat(f(:), IsEqualTo([f1(:); f2(:)], ...
      'Within', AbsoluteTolerance(tol)), ...
      'Forces do not match');
end
