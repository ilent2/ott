function tests = testSphere
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruction(testCase)

  posStiffness = 1.0;
  rotStiffness = 1.0;
  a = ott.optics.harmonic.OneAxis('posStiffness', posStiffness, ...
      'rotStiffness', rotStiffness);

  testCase.verifyEqual(a.rotStiffness, rotStiffness, 'rot Stiffness not set');
  testCase.verifyEqual(a.posStiffness, posStiffness, 'pos Stiffness not set');

end

