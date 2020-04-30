function tests = testRotateHelper
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testNargoutCheck(testCase)

  beam = ott.optics.vswf.bsc.Bsc();
  
  testCase.verifyWarning(@() beam.rotateX(5), ...
    'ott:utils:RotateHelper:nargout_check');
end
