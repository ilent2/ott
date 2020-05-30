function tests = testMaskedNearfield3d
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../../');
end

function testConstructor(testCase)

  masked_beam = ott.beam.abstract.PlaneWave();
  mask = @(x) true(size(x, 2));

  beam = ott.beam.abstract.MaskedNearfield3d(mask, masked_beam);

  testCase.verifyEqual(beam.masked_beam, masked_beam, 'beam');
  testCase.verifyEqual(beam.mask, mask, 'mask');

end

