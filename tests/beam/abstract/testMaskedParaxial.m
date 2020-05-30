function tests = testMaskedParaxial
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../../');
end

function testConstructor(testCase)

  masked_beam = ott.beam.abstract.UniformFarfield([1;0]);
  mask = @(x) true(size(x, 2));
  
  beam = ott.beam.abstract.MaskedParaxial(mask, masked_beam);
  
  testCase.verifyEqual(beam.masked_beam, masked_beam, 'beam');
  testCase.verifyEqual(beam.mask, mask, 'mask');
  testCase.verifyEqual(beam.mapping, 'tan', 'mapping');
  testCase.verifyEqual(beam.hemisphere, 'pos');

end

function testTopHat(testCase)

  background = ott.beam.abstract.UniformFarfield([1;0]);
  beam = ott.beam.abstract.MaskedParaxial.TopHat(pi/2, background);
  testCase.verifyEqual(beam.masked_beam, background, 'beam');
  
end

function testAnnular(testCase)

  background = ott.beam.abstract.UniformFarfield([1;0]);
  beam = ott.beam.abstract.MaskedParaxial.Annular([0, pi/2], background);
  testCase.verifyEqual(beam.masked_beam, background, 'beam');
  
end

