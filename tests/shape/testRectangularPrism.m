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
  testCase.verifyEqual(shape.maxRadius, sqrt(sum((widths./2).^2)), 'maxR');
  testCase.verifyEqual(shape.starShaped, true, 'star');
  testCase.verifyEqual(shape.zRotSymmetry, 2, 'zrot');
  testCase.verifyEqual(shape.xySymmetry, true, 'xy');
  testCase.verifyEqual(shape.volume, prod(widths), 'volume');
  testCase.verifyEqual(shape.boundingBox, ...
      [-1,1;-1,1;-1,1].*[0.5;1;1.5], 'bb');

  % Scale
  shape = shape * 2;
  testCase.verifyEqual(shape.widths, 2*widths, 'scaled');
  
  % Cube zRotSymmetry
  widths = [1;1;10];
  shape = ott.shape.RectangularPrism(widths);
  testCase.verifyEqual(shape.zRotSymmetry, 4, 'zrot4');

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

function testGetFaces(testCase)

  shape = ott.shape.RectangularPrism([1,1,1]);
  cube = ott.shape.Cube(shape);
  
  testCase.verifyEqual(shape.faces, cube.faces, 'cube should match RP');

end

function testMethods(testCase)

  shape = ott.shape.RectangularPrism([1,1,1]);
  
  testCase.verifyEqual(shape.normalsXyz([0.5;0;0]), [1;0;0], 'normals');
  testCase.verifyEqual(shape.insideXyz([0;0;0]), true, 'inside');
  testCase.verifyEqual(shape.intersect([0;0;0], [2;0;0]), [0.5;0;0], 'intersect');

end
