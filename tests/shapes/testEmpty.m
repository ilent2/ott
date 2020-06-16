function tests = testEmpty
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  shape = ott.shapes.Empty();

end

