function tests = testQuaternionRmatrix
  % Tests for conversion between quaternions and rotation matrices
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testForwardBack(testCase)

  rmatrix = ott.utils.roty(-0.1022);
  quat = ott.utils.rmatrix2quaternion(rmatrix);
  
  testCase.verifyEqual(size(quat), [4, 1], ...
    'Quaterinion has wrong size');
  testCase.verifyEqual(vecnorm(quat), 1, ...
    'Vecnorm of quaternion not satisfied');
  
  rmatrix2 = ott.utils.quaternion2rmatrix(quat);
  
  testCase.verifyEqual(rmatrix2, rmatrix, ...
    'AbsTol', 1e-15, ...
    'Transform back to rmatrix failed');

end

