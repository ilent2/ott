function tests = testVariable
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testFromShape(testCase)

  shape = ott.shape.Sphere(1e-6);
  ri = 1.2;
  particle = ott.particle.Variable.FromShape(...
      'initial_relative_index', ri, 'initial_shape', shape);

  % Compare to Fixed
  particleF = ott.particle.Fixed.FromShape(shape, ...
      'index_relative', ri);

  testCase.verifyEqual(particle.drag, particleF.drag, 'drag');
  testCase.verifyEqual(particle.tmatrix, particleF.tmatrix, 'tmatrix');
  testCase.verifyEqual(particle.tinternal, particleF.tinternal, 'int');
  testCase.verifyEqual(particle.shape, particleF.shape, 'shape');

end

function testSphere(testCase)

  shape = ott.shape.Sphere(1e-6);
  ri = 1.2;
  particle = ott.particle.Variable.Sphere(...
      'initial_relative_index', ri, 'initial_shape', shape);

  % Compare to Fixed
  particleF = ott.particle.Fixed.FromShape(shape, ...
      'index_relative', ri);

  testCase.verifyEqual(particle.drag, particleF.drag, 'drag');
  testCase.verifyEqual(particle.tmatrix, particleF.tmatrix, 'tmatrix');
  testCase.verifyEqual(particle.tinternal, particleF.tinternal, 'int');
  testCase.verifyEqual(particle.shape, particleF.shape, 'shape');

end

function testStarShaped(testCase)

  shape = ott.shape.Sphere(1e-6);
  ri = 1.2;
  particle = ott.particle.Variable.StarShaped(...
      'initial_relative_index', ri, 'initial_shape', shape);

end

function testSetErrors(testCase)

  shape = ott.shape.Sphere(1e-6);
  ri = 1.2;
  particle = ott.particle.Variable.FromShape(...
      'initial_relative_index', ri, 'initial_shape', shape);

  testCase.verifyError(@testSetDrag, 'ott:particle:variable:set_drag');
  testCase.verifyError(@testSetTmatrix, 'ott:particle:variable:set_tmatrix');
  testCase.verifyError(@testSetTinernal, 'ott:particle:variable:set_tinternal');
  
  function testSetDrag()
    particle.drag = 5;
  end
  function testSetTmatrix()
    particle.tmatrix = 5;
  end
  function testSetTinernal()
    particle.tinternal = 5;
  end
end

function testSetShapeRi(testCase)

  shape = ott.shape.Sphere(1e-6);
  ri = 1.2;
  particle = ott.particle.Variable.FromShape();
  
  particle.shape = shape;
  particle.relative_index = ri;

end