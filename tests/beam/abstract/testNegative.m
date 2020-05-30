function tests = testScattered
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  beam = ott.beam.abstract.Gaussian(1.0);
  nbeam = ott.beam.abstract.Negative(beam);

  testCase.verifyEqual(nbeam.data, beam);

end

