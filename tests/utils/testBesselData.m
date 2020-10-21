function tests = testBesselData
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  func = @ott.utils.sbesselj;
  data = ott.utils.BesselData(func);
  testCase.verifyEqual(data.func, func, 'incorrect function');
  testCase.verifyEmpty(data.un, 'un');
  testCase.verifyEmpty(data.ur, 'ur');
  testCase.verifyEmpty(data.hn, 'hn');
  testCase.verifyEmpty(data.dhn, 'dhn');

end

function testEvaluate(testCase)

  func = @ott.utils.sbesselj;
  data = ott.utils.BesselData(func);
  n = [1, 2, 3, 4];
  kr = linspace(0, 5, 5);
  [hn, dhn, data2] = data.evaluate(n, kr);

  testCase.verifySize(hn, [4, 5], 'hn size');
  testCase.verifySize(dhn, [4, 5], 'hn size');
  testCase.verifyEqual(data2.hn, hn, 'data2.hn');
  testCase.verifyEqual(data2.dhn, dhn, 'data2.dhn');
  testCase.verifyEqual(data2.un, n, 'data2.un');
  testCase.verifyEqual(data2.ur, kr, 'data2.ur');

end

