function tests = wavefrontobj
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testLoad(testCase)

  fn = tempname();

  % Generate a cube and write it to a file
  cube = ott.shapes.Cube(1.0);
  cube.writeWavefrontObj(fn);

  % Load the shape
  shape = ott.shapes.WavefrontObj(fn);

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.RelativeTolerance;
  tol = 1.0e-3;

  testCase.verifyThat(shape.maxRadius, IsEqualTo(cube.maxRadius, ...
      'Within', RelativeTolerance(tol)), ...
      'Incorrect radius after writing/loading OBJ file');

  testCase.verifyThat(shape.volume, IsEqualTo(cube.volume, ...
      'Within', RelativeTolerance(tol)), ...
      'Incorrect volume after writing/loading OBJ file');

  delete(fn);

end
