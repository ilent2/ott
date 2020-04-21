function tests = testSphere
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruction(testCase)

  posStiffness = zeros(3, 3);
  rotStiffness = zeros(3, 3);
  a = ott.optics.harmonic.TwoAxis('posStiffness', posStiffness, ...
      'rotStiffness1', rotStiffness, 'rotStiffness2', rotStiffness);

  testCase.verifyEqual(a.posStiffness, posStiffness, 'pos Stiffness not set');
  testCase.verifyEqual(a.rotStiffness1, rotStiffness, 'rot1 Stiffness not set');
  testCase.verifyEqual(a.rotStiffness2, posStiffness, 'rot2 Stiffness not set');

end

