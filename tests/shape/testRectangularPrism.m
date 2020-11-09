function tests = testRectangularPrism
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  widths = [1; 2; 3];
  shape = ott.shape.RectangularPrism(widths);

  testCase.verifyEqual(shape.widths, widths, 'widths');
  testCase.verifyEqual(shape.starShaped, true, 'star');
  testCase.verifyEqual(shape.zRotSymmetry, 2, 'zrot');
  testCase.verifyEqual(shape.xySymmetry, true, 'xy');
  testCase.verifyEqual(shape.volume, prod(widths), 'volume');
  testCase.verifyEqual(shape.boundingBox, ...
      [-1,1;-1,1;-1,1].*[0.5;1;1.5], 'bb');

  % Scale
  shape = shape * 2;
  testCase.verifyEqual(shape.widths, 2*widths, 'scaled');

end

function testSetErrors(testCase)

  widths = [1; 2; 3];
  shape = ott.shape.RectangularPrism(widths);

  testCase.verifyError(@setRadius, 'ott:shape:RectangularPrism:set_maxradius');
  testCase.verifyError(@setVolume, 'ott:shape:RectangularPrism:set_volume');

  function setRadius
    shape.maxRadius = 1;
  end

  function setVolume
    shape.volume = 1;
  end
end

