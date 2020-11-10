function tests = testPmParaxial
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstructBeam(testCase)

  truncA = pi/4;
  profile = ott.beam.Gaussian('truncation_angle', truncA);
  profile.data = profile.data.setNmax(10, 'RelTol', []);
  beam = ott.beam.PmParaxial.BeamProfile(profile, ...
    'Nmax', profile.data.Nmax, 'direction', 'neg', ...
    'truncation_angle', truncA);

  testCase.verifyEqual(beam.data.power, profile.data.power, ...
      'RelTol', 0.15, 'power');

end

function testConstructInterp1d(testCase)
  % Construct using a set of points

  x = linspace(-1, 1, 20);
  y = linspace(-1, 1, 20);
  [X, Y] = ndgrid(x, y);
  amp = 2 - (X.^2 + Y.^2);

  polbasis = 'cartesian';
  polfield = [1, 0];
  Nmax = 1;
  beam = ott.beam.PmParaxial.InterpProfile(X, Y, amp, ...
      'polbasis', polbasis, 'polfield', polfield, ...
      'Nmax', Nmax);
    
  testCase.verifyEqual(beam.data.Nmax, Nmax, 'Nmax');
  testCase.verifyEqual(beam.polfield, polfield, 'polfield');
  testCase.verifyEqual(beam.polbasis, polbasis, 'polbasis');

end

function testConstructInterp2d(testCase)
  % Construct using a set of points

  x = linspace(-1, 1, 20);
  y = linspace(-1, 1, 20);
  [X, Y] = ndgrid(x, y);
  amp = 2 - (X.^2 + Y.^2);
  amp = permute(cat(3, amp, amp), [3, 1, 2]);

  polbasis = 'polar';
  polfield = [1, 0];
  mapping = 'tan';
  Nmax = 1;
  beam = ott.beam.PmParaxial.InterpProfile(X, Y, amp, ...
      'polbasis', polbasis, 'polfield', polfield, ...
      'Nmax', Nmax, 'mapping', mapping);
    
  testCase.verifyEqual(beam.data.Nmax, Nmax, 'Nmax');
  testCase.verifyEqual(beam.polfield, polfield, 'polfield');
  testCase.verifyEqual(beam.polbasis, polbasis, 'polbasis');

end

