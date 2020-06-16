function tests = testRectangularPrism
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  widths = [1; 2; 3];
  shape = ott.shapes.RectangularPrism(widths);

  testCase.verifyEqual(shape.widths, widths, 'widths');

end

