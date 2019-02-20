function tests = testCol3to1
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testBasic(testCase)

  import matlab.unittest.constraints.IsEqualTo;

  ret = ott.utils.col3to1([1, 2, 3; 4, 5, 6; 7, 8, 9]);
  testCase.verifyThat(ret, IsEqualTo((1:9).'), ...
      'Failed col3to1');

  ret = ott.utils.col1to3([1; 2; 3; 4; 5; 6; 7; 8; 9]);
  testCase.verifyThat(ret, IsEqualTo([1, 2, 3; 4, 5, 6; 7, 8, 9]), ...
      'Failed col1to3');

end

