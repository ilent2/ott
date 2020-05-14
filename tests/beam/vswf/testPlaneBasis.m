function tests = testPlaneBasis
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructOptional(testCase)

  direction = randn(3, 2);
  polarisation = randn(3, 2);
  Nmax = 1;

  beam = ott.beam.vswf.PlaneBasis(Nmax, 'direction', direction, ...
      'polarisation', polarisation);

  testCase.verifyEqual(numel(beam), 2, 'numel');
  testCase.verifyEqual(beam.direction, direction, 'direction');
  testCase.verifyEqual(beam.polarisation, polarisation, 'polarisation');
  testCase.verifyEqual(beam.Nmax, Nmax, 'nmax');

end

function testTranslateZ(testCase)
  % Verify translations simply shift the phases of the plane waves

  assert(false);

end

