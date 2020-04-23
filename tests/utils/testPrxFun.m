function tests = testPrxFun
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testRotation(testCase)

  rotFun = @(~, pos, ~, rot) rot * pos;
  
  orot = ott.utils.rotz(90);
  opos = rand(3, 5);
  
  out = ott.utils.prxfun(rotFun, 3, 'position', opos, 'rotation', orot);
  
  testCase.verifyEqual(out, orot * opos, ...
    'AbsTol', 1e-15, 'Position of arguments doesnt work');
end
