function tests = testLgParaxialBasis
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

% See testConsistency for default constructor

function testConstructOptional(testCase)

  lmodes = [0, 1, 2];
  amodes = [0, 2, 4];
  paraxial_waist = 1.0;

  beam = ott.beam.vswf.LgParaxialBasis(paraxial_waist, lmodes, amodes);

  testCase.verifyEqual(beam.paraxial_waist, paraxial_waist, 'waist');
  testCase.verifyEqual(beam.lmodes, lmodes, 'lmodes');
  testCase.verifyEqual(beam.amodes, amodes, 'amodes');
end

function testAngularScaling(testCase)

  beam = ott.optics.vswf.bsc.PmGauss('angular_scaling', 'sintheta');
  testCase.verifyEqual(beam.angular_scaling, 'sintheta', ...
    'sintheta Angular scaling not set correctly');

  beam = ott.optics.vswf.bsc.PmGauss('angular_scaling', 'tantheta');
  testCase.verifyEqual(beam.angular_scaling, 'tantheta', ...
    'tantheta Angular scaling not set correctly');

end

