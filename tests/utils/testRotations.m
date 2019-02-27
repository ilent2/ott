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

end

function testRoty(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-4;

  roty45 = [0.7071         0    0.7071
                 0    1.0000         0
           -0.7071         0    0.7071];

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

end

function testRotx(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-4;

  rotx45 = [1.0000         0         0
                 0    0.7071   -0.7071
                 0    0.7071    0.7071];

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

end

