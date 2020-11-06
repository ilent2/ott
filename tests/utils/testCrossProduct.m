function tests = testCrossProduct
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testWithCrossProduct(testCase)

  a = [1;2;3];
  b = [4;5;6];
  
  target = cross(a, b);
  
  mat = ott.utils.crossProductMatrix(a);
  trial = mat * b;
  
  testCase.verifyEqual(trial, target, 'RelTol', 1e-15);

end

