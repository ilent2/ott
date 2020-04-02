function tests = testVswfCart
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testFiniteAtZero(testCase)

  import matlab.unittest.constraints.IsFinite

  rtp = [0, 0, 0];

  [A, B] = ott.utils.vswfcart(1, 0, rtp(1), rtp(2), rtp(3), 'regular');
  testCase.verifyThat(A, IsFinite, 'A not finite');
  testCase.verifyThat(B, IsFinite, 'B not finite');

end

