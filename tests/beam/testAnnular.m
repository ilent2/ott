function tests = testAnnular
  tests = functiontests(localfunctions)
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  theta = [2*pi/8, 3*pi/8];
  polbasis = 'polar';
  polfield = [1, 0];
  profile = 'uniform';
  Nmax = 2;
  beam = ott.beam.Annular(theta, profile, ...
      'polbasis', polbasis, 'polfield', polfield, ...
      'initial_Nmax', Nmax);

  testCase.verifyEqual(beam.theta, theta, 'theta');
  testCase.verifyEqual(beam.profile, profile, 'profile');
  testCase.verifyEqual(beam.polbasis, polbasis, 'polbasis');
  testCase.verifyEqual(beam.polfield, polfield, 'polfield');
  testCase.verifyEqual(beam.data.Nmax, Nmax, 'initial_Nmax');

end

function testConstructGaussian(testCase)
  % Construct a Gaussian masked by a annular

  theta = [0, pi];   % full range, no masking
  profile = ott.beam.Gaussian();
  beam = ott.beam.Annular(theta, profile, 'initial_Nmax', profile.Nmax);

  testCase.verifyEqual(beam.profile, profile, 'profile');
  testCase.verifyEqual(beam.data.a, profile.a, 'beam a');
  testCase.verifyEqual(beam.data.b, profile.b, 'beam b');

end

function testConstructInterp(testtCase)
  % Construct using a set of points

  theta = linspace(0, pi/2);
  amp = linspace(0, 1);
  polbasis = 'polar';
  polfield = [1, 0];
  Nmax = 2;
  beam = ott.beam.Annular.InterpProfile(theta, amp, ...
      'polbasis', polbasis, 'polfield', polfield, ...
      'initial_Nmax', Nmax);

end

