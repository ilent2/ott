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

function testContains(testCase)

  beam1 = ott.beam.abstract.Gaussian(1.0);
  array_beam = ott.beam.Array('incoherent', {beam1});
  beam = ott.beam.abstract.Negative(array_beam);

  testCase.verifyEqual(beam.contains('incoherent'), true, '1incoherent');
  testCase.verifyEqual(beam.contains('coherent'), false, '1coherent');
  testCase.verifyEqual(beam.contains('array'), false, '1array');

end

