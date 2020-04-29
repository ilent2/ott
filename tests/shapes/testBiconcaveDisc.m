function tests = testBiconcaveDisc
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  % Test defaults
  a = [0.0518, 2.0026, -4.491];
  radius = 1.0;
  shape = ott.shapes.BiconcaveDisc();
  
  testCase.verifyEqual(shape.radius, radius, 'default rad');
  testCase.verifyEqual(shape.coefficients, a, 'default a');

  % Test args
  a = 2*[0.0518, 2.0026, -4.491];
  radius = 2.0;
  shape = ott.shapes.BiconcaveDisc(radius, a);
  
  testCase.verifyEqual(shape.radius, radius, 'coeff rad');
  testCase.verifyEqual(shape.coefficients, a, 'coeff a');
  
  testCase.verifyGreaterThan(shape.maxRadius, radius, 'maxRadius');
  testCase.verifyEqual(shape.maxXyRadius, radius, 'maxXyRadius');
  testCase.verifyEqual(size(shape.volume), [1, 1], 'vol')

end

function testRadialProfile(testCase)

  radius = 1.0;
  r = linspace(0, 2*radius, 100);
  shape = ott.shapes.BiconcaveDisc(radius);
  z = shape.extrudeProfile(r);
  
  testCase.verifyEqual(size(z), size(r), 'size change');
  testCase.verifyEqual(z(r>radius), nan*z(r>radius), 'nans missing');

end

function testSurf(testCase)

  shape = ott.shapes.BiconcaveDisc();
  
  h = figure();
  shape.surf('show_normals', true);
  axis equal;
  close(h);
  
end
