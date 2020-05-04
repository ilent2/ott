function tests = testSphere
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruction(testCase)

  stiffness = 1.0;
  a = ott.optics.harmonic.Sphere(stiffness);

  testCase.verifyEqual(a.stiffness, stiffness, 'Stiffness not set');

end

