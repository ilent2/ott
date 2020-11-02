function tests = testPlaneWave
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  polarisation = [1, 0];
  Nmax = 2;
  beam = ott.beam.PlaneWave(polarisation, 'initial_Nmax', Nmax);

  testCase.verifyEqual(beam.data.Nmax, Nmax, 'Nmax');
  testCase.verifyEqual(beam.polarisation, polarisation, 'pol');
  testCase.verifyEqual(beam.data.direction, [0;0;1], 'direction');

end

