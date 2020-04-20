function tests = testRotatePolarizability
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testFunction(testCase)

  import matlab.unittest.constraints.IsEqualTo;

  alpha = eye(3);
  rotation = eye(3);

  A1 = ott.utils.rotate_polarizability(alpha, rotation);
  testCase.verifyThat(size(A1), IsEqualTo(size(alpha)), ...
      'Incorrect size for matrix A1');

  A3 = ott.utils.rotate_polarizability(alpha, 'sphdirection', [0.1; 0.5]);
  testCase.verifyThat(size(A3), IsEqualTo(size(alpha)), ...
      'Incorrect size for matrix A3');

end

function testDirection(testCase)

  alpha = diag([1, 1, 3]);
  result = diag([1, 3, 1]);
  direction = [0; 1; 0];
  A = ott.utils.rotate_polarizability(alpha, 'dir', direction);
  testCase.verifyEqual(A, result, 'AbsTol', 1e-6, ...
    'Rotation to y-axis failed');

  alpha = diag([1, 1, 3]);
  result = diag([3, 1, 1]);
  direction = [1; 0; 0];
  A = ott.utils.rotate_polarizability(alpha, 'dir', direction);
  testCase.verifyEqual(A, result, 'AbsTol', 1e-6, ...
    'Rotation to x-axis failed');

  alpha = diag([1, 1, 3]);
  result = diag([1, 1, 3]);
  direction = [0; 0; 1];
  A = ott.utils.rotate_polarizability(alpha, 'dir', direction);
  testCase.verifyEqual(A, result, 'AbsTol', 1e-6, ...
    'Rotation to z-axis failed');
  
end

function testDirectionArray(testCase)

  alpha = diag([1, 1, 3]);
  direction = [1, 0, 0; 0, 1, 0; 0, 0, 1].';
  result = [diag([3, 1, 1]), diag([1, 3, 1]), diag([1, 1, 3])];
  A = ott.utils.rotate_polarizability(alpha, 'dir', direction);
  testCase.verifyEqual(A, result, 'AbsTol', 1e-6, ...
    'Multiple directions with single alpha failed');
  

end
