function tests = testFixed
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testFromShape(testCase)

  shape = ott.shape.Sphere(1e-6);
  index_relative = 1.2;
  particle = ott.particle.Fixed.FromShape(shape, ...
      'index_relative', index_relative, 'internal', true);

  testCase.verifyEqual(particle.drag.viscosity, 8.9e-4, 'viscosity');
  testCase.verifyEqual(particle.tmatrix.relative_index, index_relative, 'ri');
  testCase.verifyEqual(particle.tinternal.relative_index, index_relative, 'tinternal');
  testCase.verifyEqual(particle.shape, shape, 'shape');

end

