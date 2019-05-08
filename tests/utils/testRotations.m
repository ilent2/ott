function tests = testRotations
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testRotz(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-4;

  rotz45 = [0.7071   -0.7071         0
            0.7071    0.7071         0
                 0         0    1.0000];
               
  rotz90= [0, -1, 0; 1, 0, 0; 0, 0, 1];
  rotz180= [-1, 0, 0; 0, -1, 0; 0, 0, 1];

  val = ott.utils.rotz(45);
  testCase.verifyThat(val, IsEqualTo(rotz45, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect rotz values');

  val = ott.utils.rotz([1, 2, 3; 4, 5, 6]);
  testCase.verifyThat(size(val), IsEqualTo([3, 3*6]), ...
      'Incorrect size for rotz matrix input');

  val = ott.utils.rotz([1, 2, 3; 4, 5, 6], 'usecell', true);
  testCase.verifyThat(size(val), IsEqualTo([2, 3]), ...
      'Incorrect size for rotz matrix input cell output');

  val = ott.utils.rotz([45, 90; 0, 180]);
  res = [rotz45, eye(3), rotz90, rotz180];
  testCase.verifyThat(val, IsEqualTo(res, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect output for rotz matrix input');

  val = ott.utils.rotz([45, 90; 0, 180], 'usecell', true);
  res = {rotz45, rotz90; eye(3), rotz180};
  testCase.verifyThat(val, IsEqualTo(res, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect output for rotz cell input');

end

function testRoty(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-4;

  roty45 = [0.7071         0    0.7071
                 0    1.0000         0
           -0.7071         0    0.7071];
               
  roty90= [0, 0, 1; 0, 1, 0; -1, 0, 0];
  roty180= [-1, 0, 0; 0, 1, 0; 0, 0, -1];

  val = ott.utils.roty(45);
  testCase.verifyThat(val, IsEqualTo(roty45, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect roty values');

  val = ott.utils.roty([1, 2, 3; 4, 5, 6]);
  testCase.verifyThat(size(val), IsEqualTo([3, 3*6]), ...
      'Incorrect size for roty matrix input');

  val = ott.utils.roty([1, 2, 3; 4, 5, 6], 'usecell', true);
  testCase.verifyThat(size(val), IsEqualTo([2, 3]), ...
      'Incorrect size for roty matrix input cell output');

  val = ott.utils.roty([45, 90; 0, 180]);
  res = [roty45, eye(3), roty90, roty180];
  testCase.verifyThat(val, IsEqualTo(res, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect output for roty matrix input');

  val = ott.utils.roty([45, 90; 0, 180], 'usecell', true);
  res = {roty45, roty90; eye(3), roty180};
  testCase.verifyThat(val, IsEqualTo(res, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect output for roty cell input');

end

function testRotx(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-4;

  rotx45 = [1.0000         0         0
                 0    0.7071   -0.7071
                 0    0.7071    0.7071];
               
  rotx90= [1, 0, 0; 0, 0, -1; 0, 1, 0];
  rotx180= [1, 0, 0; 0, -1, 0; 0, 0, -1];

  val = ott.utils.rotx(45);
  testCase.verifyThat(val, IsEqualTo(rotx45, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect rotx values');

  val = ott.utils.rotx([1, 2, 3; 4, 5, 6]);
  testCase.verifyThat(size(val), IsEqualTo([3, 3*6]), ...
      'Incorrect size for rotx matrix input');

  val = ott.utils.rotx([1, 2, 3; 4, 5, 6], 'usecell', true);
  testCase.verifyThat(size(val), IsEqualTo([2, 3]), ...
      'Incorrect size for rotx matrix input cell output');

  val = ott.utils.rotx([45, 90; 0, 180]);
  res = [rotx45, eye(3), rotx90, rotx180];
  testCase.verifyThat(val, IsEqualTo(res, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect output for rotx matrix input');

  val = ott.utils.rotx([45, 90; 0, 180], 'usecell', true);
  res = {rotx45, rotx90; eye(3), rotx180};
  testCase.verifyThat(val, IsEqualTo(res, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect output for rotx cell input');

end

