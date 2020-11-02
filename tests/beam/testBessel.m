function tests = testBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  theta = pi/4;
  polfield = [1, 0];
  lmode = 1;
  polbasis = 'polar';
  Nmax = 2;
  beam = ott.beam.Bessel(theta, polfield, lmode, 'polbasis', polbasis, ...
      'initial_Nmax', Nmax);

  testCase.verifyEqual(beam.polfield, polfield, 'polfield');
  testCase.verifyEqual(beam.polbasis, polfield, 'polbasis');
  testCase.verifyEqual(beam.data.Nmax, Nmax, 'Nmax');
  testCase.verifyEqual(beam.lmode, lmode, 'lmode');
  testCase.verifyEqual(beam.theta, theta, 'theta');

end

