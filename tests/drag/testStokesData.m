function tests = testStokes
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  data = rand(6, 6);

  drag = ott.drag.StokesData(data);
  testCase.verifyEqual(drag.forward, data, 'forward');

end

