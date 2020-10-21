function tests = testMatchsize
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testSimple(testCase)

  a = zeros(1, 5);
  b = zeros(5, 6);
  [A, B] = ott.utils.matchsize(a, b);

  testCase.verifySize(A, [5, 5], 'incorrect A size');
  testCase.verifySize(B, [5, 6], 'incorrect B size');

end

