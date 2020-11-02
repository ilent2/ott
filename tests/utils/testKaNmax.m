function tests = testKaNmax
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testZeroNmax(testCase)

  testCase.verifyEqual(ott.utils.ka2nmax(0), 0, 'ka2nmax');
  testCase.verifyEqual(ott.utils.nmax2ka(0), 0, 'nmax2ka');

end

