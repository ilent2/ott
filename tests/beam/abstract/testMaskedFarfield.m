function tests = testMaskedFarfield
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../../');
end

function testConstructor(testCase)

  masked_beam = ott.beam.abstract.UniformFarfield([1;0]);
  mask = @(x) true(size(x, 2));
  
  beam = ott.beam.abstract.MaskedFarfield(mask, masked_beam);
  
  testCase.verifyEqual(beam.masked_beam, masked_beam, 'beam');
  testCase.verifyEqual(beam.mask, mask, 'mask');

end

function testContains(testCase)

  beam1 = ott.beam.abstract.Gaussian(1.0);
  masked_beam = ott.beam.Array('incoherent', {beam1});
  mask = @(x) true(size(x, 2));
  beam = ott.beam.abstract.MaskedFarfield(mask, masked_beam);

  testCase.verifyEqual(beam.contains('incoherent'), true, '1incoherent');
  testCase.verifyEqual(beam.contains('coherent'), false, '1coherent');
  testCase.verifyEqual(beam.contains('array'), false, '1array');

  beam1 = ott.beam.abstract.Gaussian(1.0);
  masked_beam = ott.beam.Array('coherent', {beam1});
  mask = @(x) true(size(x, 2));
  beam = ott.beam.abstract.MaskedFarfield(mask, masked_beam);

  testCase.verifyEqual(beam.contains('incoherent'), false, '2incoherent');
  testCase.verifyEqual(beam.contains('coherent'), true, '2coherent');
  testCase.verifyEqual(beam.contains('array'), false, '2array');

end

function testTopHat(testCase)

  background = ott.beam.abstract.UniformFarfield([1;0]);
  beam = ott.beam.abstract.MaskedFarfield.TopHat(pi/2, background);
  testCase.verifyEqual(beam.masked_beam, background, 'beam');
  
end

function testAnnular(testCase)

  background = ott.beam.abstract.UniformFarfield([1;0]);
  beam = ott.beam.abstract.MaskedFarfield.Annular([0, pi/2], background);
  testCase.verifyEqual(beam.masked_beam, background, 'beam');
  
end

