function tests = testWignerRotationMatrix
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testEye(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  % Test Nmax = 1
  val = ott.utils.wigner_rotation_matrix(1, eye(3));
  testCase.verifyThat(full(val), IsEqualTo(eye(3), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect results for Nmax=1 R=eye(3)');

  % Test Nmax = 2
  Nmax = 2;
  val = ott.utils.wigner_rotation_matrix(Nmax, eye(3));
  numVals = ott.utils.combined_index(Nmax, Nmax);
  testCase.verifyThat(val, IsEqualTo(speye(numVals), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect results for Nmax=2 R=eye(3)');

  % Test Nmax = 20
  Nmax = 20;
  val = ott.utils.wigner_rotation_matrix(Nmax, eye(3));
  numVals = ott.utils.combined_index(Nmax, Nmax);
  testCase.verifyThat(val, IsEqualTo(speye(numVals), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect results for Nmax=200 R=eye(3)');

end

function testArbitraryValue1(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-4;

  % Results for Nmax=3, R = rotx(32)*rotz(18))
  n1 = [0.8788 - 0.2855i   0.1158 + 0.3564i  -0.0723 + 0.0235i;
        0.0000 + 0.3747i   0.8480 + 0.0000i   0.0000 + 0.3747i;
       -0.0723 - 0.0235i  -0.1158 + 0.3564i   0.8788 + 0.2855i];
  n2 = [0.6908 - 0.5019i   0.2878 + 0.3961i  -0.1391 + 0.1011i  -0.0237 - 0.0326i   0.0047 - 0.0034i;
        0.1513 + 0.4657i   0.6117 - 0.1988i   0.1701 + 0.5235i  -0.1948 + 0.0633i  -0.0124 - 0.0383i;
       -0.1720 + 0.0000i   0.0000 + 0.5504i   0.5788 + 0.0000i   0.0000 + 0.5504i  -0.1720 + 0.0000i;
        0.0124 - 0.0383i  -0.1948 - 0.0633i  -0.1701 + 0.5235i   0.6117 + 0.1988i  -0.1513 + 0.4657i;
        0.0047 + 0.0034i   0.0237 - 0.0326i  -0.1391 - 0.1011i  -0.2878 + 0.3961i   0.6908 + 0.5019i];
  n3 = [0.4637 - 0.6383i   0.4483 + 0.3257i  -0.1477 + 0.2033i  -0.0673 - 0.0489i   0.0121 - 0.0167i   0.0030 + 0.0022i  -0.0003 + 0.0004i;
        0.3257 + 0.4483i   0.3759 - 0.2731i   0.3513 + 0.4836i  -0.2638 + 0.1917i  -0.0663 - 0.0913i   0.0212 - 0.0154i   0.0022 + 0.0030i;
       -0.2389 + 0.0776i   0.1847 + 0.5685i   0.2872 - 0.0933i   0.1841 + 0.5665i  -0.3300 + 0.1072i  -0.0349 - 0.1073i   0.0196 - 0.0064i;
        0.0000 - 0.0832i  -0.3261 + 0.0000i   0.0000 + 0.5957i   0.2527 + 0.0000i   0.0000 + 0.5957i  -0.3261 + 0.0000i   0.0000 - 0.0832i;
        0.0196 + 0.0064i   0.0349 - 0.1073i  -0.3300 - 0.1072i  -0.1841 + 0.5665i   0.2872 + 0.0933i  -0.1847 + 0.5685i  -0.2389 - 0.0776i;
       -0.0022 + 0.0030i   0.0212 + 0.0154i   0.0663 - 0.0913i  -0.2638 - 0.1917i  -0.3513 + 0.4836i   0.3759 + 0.2731i  -0.3257 + 0.4483i;
       -0.0003 - 0.0004i  -0.0030 + 0.0022i   0.0121 + 0.0167i   0.0673 - 0.0489i  -0.1477 - 0.2033i  -0.4483 + 0.3257i   0.4637 + 0.6383i];

  Nmax = 3;
  val = ott.utils.wigner_rotation_matrix(Nmax, ott.utils.rotx(32)*ott.utils.rotz(18));

  testCase.verifyThat(full(val(1:3, 1:3)), IsEqualTo(n1, ...
      'Within', AbsoluteTolerance(tol)), ...
      'n1 values disagree');
  testCase.verifyThat(full(val(4:8, 4:8)), IsEqualTo(n2, ...
      'Within', AbsoluteTolerance(tol)), ...
      'n2 values disagree');
  testCase.verifyThat(full(val(9:end, 9:end)), IsEqualTo(n3, ...
      'Within', AbsoluteTolerance(tol)), ...
      'n3 values disagree');

end

