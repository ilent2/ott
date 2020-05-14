function tests = testPointmatch
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testInversionCoefficients(testCase)

  % We want to be able to keep inverse coefficients for faster
  % repeated calculations.  We had tests for this in 671899bbee931

  assert(false);

end

