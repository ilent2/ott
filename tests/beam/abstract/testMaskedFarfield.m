function tests = testBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  masked_beam = ott.beam.abstract.UniformFarfield([1;0]);
  mask = @(x) true(size(x, 2));
  
  beam = ott.beam.abstract.MaskedFarfield(mask, masked_beam);
  
  testCase.verifyEqual(beam.beam, masked_beam, 'beam');
  testCase.verifyEqual(beam.mask, mask, 'mask');

end

function testTopHat(testCase)

  background = ott.beam.abstract.UniformFarfield([1;0]);
  beam = ott.beam.abstract.MaskedFarfield.TopHat(pi/2, background);
  testCase.verifyEqual(beam.beam, background, 'beam');
  
end

function testAnnular(testCase)

  background = ott.beam.abstract.UniformFarfield([1;0]);
  beam = ott.beam.abstract.MaskedFarfield.Annular([0, pi/2], background);
  testCase.verifyEqual(beam.beam, background, 'beam');
  
end

