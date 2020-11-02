function tests = testBscBeam
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testGetData(testCase)

  beam(1) = ott.beam.Gaussian();
  beam(2) = ott.beam.LaguerreGaussian();

  Nmax = 5;
  bsc = beam.getData(Nmax);

  testCase.assertSize(bsc, size(beam), 'size');
  testCase.verifyEqual(bsc(1), beam(1).data, '1');
  testCase.verifyEqual(bsc(2), beam(2).data, '2');

end

