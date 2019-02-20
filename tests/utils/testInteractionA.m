function tests = interaction_A
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testFunction(testCase)

  import matlab.unittest.constraints.IsEqualTo;

  k = 2*pi;
  voxels = [1, 2, 3; 4, 5, 6];
  alpha = 2.0;

  A1 = ott.utils.interaction_A(k, voxels, alpha);
  testCase.verifyThat(size(A1), IsEqualTo([1,1].*size(voxels, 1)*3), ...
      'Incorrect size for matrix A1');

  A2 = ott.utils.interaction_A(k, voxels, 'inv_alpha', 1.0/alpha);
  testCase.verifyThat(size(A1), IsEqualTo([1,1].*size(voxels, 1)*3), ...
      'Incorrect size for matrix A2');

  testCase.verifyThat(A1, IsEqualTo(A2), ...
      'alpha and inv_alpha isotropic do not match');

  alpha = eye(3)*2;
  A3 = ott.utils.interaction_A(k, voxels, alpha);
  testCase.verifyThat(size(A3), IsEqualTo([1,1].*size(voxels, 1)*3), ...
      'Incorrect size for matrix A2');

  A4 = ott.utils.interaction_A(k, voxels, 'inv_alpha', inv(alpha));
  testCase.verifyThat(size(A4), IsEqualTo([1,1].*size(voxels, 1)*3), ...
      'Incorrect size for matrix A2');

end

