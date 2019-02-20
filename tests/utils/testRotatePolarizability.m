function tests = testRotatePolarizability
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testFunction(testCase)

  import matlab.unittest.constraints.IsEqualTo;

  alpha = eye(3);
  rotation = eye(3);
  direction = [1; 0; 0];

  A1 = ott.utils.rotate_polarizability(alpha, rotation);
  testCase.verifyThat(size(A1), IsEqualTo(size(alpha)), ...
      'Incorrect size for matrix A1');

  A2 = ott.utils.rotate_polarizability(alpha, 'direction', direction);
  testCase.verifyThat(size(A2), IsEqualTo(size(alpha)), ...
      'Incorrect size for matrix A2');

  A3 = ott.utils.rotate_polarizability(alpha, 'sphdirection', [0.1; 0.5]);
  testCase.verifyThat(size(A3), IsEqualTo(size(alpha)), ...
      'Incorrect size for matrix A3');

end

