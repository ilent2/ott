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
  Nmax = 20;
  beam = ott.beam.Annular(theta, ...
      'polbasis', polbasis, 'polfield', polfield, ...
      'Nmax', Nmax);

  testCase.verifyEqual(beam.theta, theta, 'theta');
  testCase.verifyEqual(beam.polbasis, polbasis, 'polbasis');
  testCase.verifyEqual(beam.polfield, polfield, 'polfield');
  testCase.verifyEqual(beam.Nmax, Nmax, 'iNmax');

end

function testConstructGaussian(testCase)
  % Construct a Gaussian masked by a annular

  theta = [0, pi];   % full range, no masking
  profile = ott.beam.Gaussian();
  beam = ott.beam.Annular.BeamProfile(theta, profile, ...
      'Nmax', profile.data.Nmax);

  bsc = ott.bsc.Bsc(beam, profile.data.Nmax);
  testCase.verifyEqual(bsc.power, profile.data.power, ...
      'RelTol', 0.2, 'power');

end

function testConstructInterp(testCase)
  % Construct using a set of points

  theta = linspace(0, pi/2);
  amp = linspace(0, 1);
  polbasis = 'polar';
  polfield = [1, 0];
  Nmax = 2;
  beam = ott.beam.Annular.InterpProfile(theta, amp, ...
      'polbasis', polbasis, 'polfield', polfield, ...
      'Nmax', Nmax);

end

