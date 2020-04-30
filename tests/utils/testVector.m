function tests = testVector
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testRotateVec(testCase)

  vec = ott.utils.Vector([0;0;1]);
  rvec = vec.rotateY(pi/2);
  
  testCase.verifyEqual(rvec.direction, [1;0;0], 'AbsTol', 1e-15, 'dir');
  testCase.verifyEqual(rvec.origin, vec.origin, 'AbsTol', 1e-15, 'origin');

end

function testRotateOrigin(testCase)

  origin = [1;0;0];
  vec = ott.utils.Vector(origin, [0;0;1]);
  rvec = vec.rotateY(pi/2, 'origin', true);
  
  testCase.verifyEqual(rvec.direction, [1;0;0], 'AbsTol', 1e-15, 'dir');
  testCase.verifyEqual(rvec.origin, [0;0;-1], 'AbsTol', 1e-15, 'origin');

end

function testScalarMul(testCase)

  origin = [1;0;0];
  vec = ott.utils.Vector(origin, [0;0;1]);
  
  svec = 2.0 * vec;
  testCase.verifyEqual(svec.direction, 2*vec.direction, 'AbsTol', 1e-15, 'dir');
  testCase.verifyEqual(svec.origin, vec.origin, 'AbsTol', 1e-15, 'origin');
  
  svec = vec * 2.0;
  testCase.verifyEqual(svec.direction, 2*vec.direction, 'AbsTol', 1e-15, 'dir');
  testCase.verifyEqual(svec.origin, vec.origin, 'AbsTol', 1e-15, 'origin');

end

