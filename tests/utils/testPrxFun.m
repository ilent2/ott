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
  
  out = ott.utils.prxfun(rotFun, 'position', opos, 'rotation', orot);
  out = cat(2, out{:});
  
  testCase.verifyEqual(out, orot * opos, ...
    'AbsTol', 1e-15, 'Position of arguments doesnt work');
end

function testMultipleOutputs(testCase)

  testData = randn(3, 5);
  testFun = @(~, pos, ~, rot) deal(pos(1), pos(2), pos(3));
  
  [out1, out2, out3] = ott.utils.prxfun(testFun, 'position', testData);
  out1 = cat(2, out1{:});
  out2 = cat(2, out2{:});
  out3 = cat(2, out3{:});
  
  testCase.verifyEqual(out1, testData(1, :), 'x');
  testCase.verifyEqual(out2, testData(2, :), 'y');
  testCase.verifyEqual(out3, testData(3, :), 'z');
  
end
