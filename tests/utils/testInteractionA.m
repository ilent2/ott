function tests = testInteractionA
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
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

function testSingleDipole(testCase)

  import matlab.unittest.constraints.IsFinite

  xyz = [0, 0, 0];
  k = 2*pi;
  
  alpha = 1.0;
  %alpha = 0.0;  % This is singular!  So is [0, 1, 1]!
  
  A1 = ott.utils.interaction_A(k, xyz, alpha);
  testCase.verifyThat(A1, IsFinite);
  
  A2 = ott.utils.interaction_A(k, xyz, [1,1,1]*alpha);
  testCase.verifyEqual(A2, A1);
  
  A3 = ott.utils.interaction_A(k, xyz, [alpha, 2, 3]);
  testCase.verifyEqual(A3(1, :), A1(1, :));
  testCase.verifyNotEqual(A3(2:3, :), A1(2:3, :));
  
end

