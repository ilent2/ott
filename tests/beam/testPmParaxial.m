function tests = testPmParaxial
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstructBeam(testCase)

  profile = ott.beam.Gaussian();
  beam = ott.beam.PmParaxial(profile);

  testCase.verifyEqual(beam.getData(10), profile.getData(10), 'data');

end

function testConstructUniform(testCase)

  theta = pi/2;
  beam = ott.beam.PmParaxial('truncation_angle', theta);
  testCase.verifyEqual(beam.truncation_angle, theta, 'theta');
end

function testConstructInterp(testCase)
  % Construct using a set of points

  x = linspace(-1, 1, 20);
  y = linspace(-1, 1, 20);
  [X, Y] = ndgrid(x, y);
  amp = X.^2 + Y.^2;

  polbasis = 'polar';
  polfield = [1, 0];
  Nmax = 1;
  beam = ott.beam.PmParaxial.InterpProfile(X, Y, amp, ...
      'polbasis', polbasis, 'polfield', polfield, ...
      'Nmax', Nmax);

end

