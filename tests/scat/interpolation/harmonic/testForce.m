function tests = testForce
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../../');
end

function testConstruction(testCase)

  stiffness = -0.5;
  position = [0;0;1] .* linspace(-5, 5, 100);
  force = stiffness .* position;

  part = ott.scat.interp.harmonic.Force(position, force);

  testCase.verifyEqual(part.stiffness, stiffness, 'AbsTol', 1e-15, 'stiffness');

end

