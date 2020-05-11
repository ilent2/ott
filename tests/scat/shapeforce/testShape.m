function tests = testShape
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruct(testCase)

  shape = ott.shapes.Sphere(1.0);
  index_relative = 1.33;
  
  particle = ott.scat.shapeforce.Shape(shape, index_relative);
  
  testCase.verifyEqual(particle.shape, shape, 'shape');
  testCase.verifyEqual(particle.index_relative, index_relative, 'ir');

end

% TODO: Test against Mie result

function testForceIndexMatched(testCase)

  beam = ott.beam.PlaneWave();
  
  shape = ott.shapes.Sphere(1.0);
  index_relative = 1.0;
  particle = ott.scat.shapeforce.Shape(shape, index_relative);
  
  force = particle.force(beam);
  testCase.verifyEqual(force, [0;0;0], 'AbsTol', 1e-15, 'index matched');
  
end

function testForceBeamArray(testCase)

  beam = ott.beam.PlaneWave('origin', randn(3, 2));
  beam.array_type = 'array';
  
  shape = ott.shapes.Sphere(1.0);
  index_relative = 1.0;
  particle = ott.scat.shapeforce.Shape(shape, index_relative);
  
  force = particle.force(beam);
  testCase.verifyEqual(force, [0;0;0].*[1,1], 'AbsTol', 1e-15, 'index matched');
  
end

function testForceBeamPosition(testCase)

  beam = ott.beam.abstract.Gaussian(1.0);
  
  shape = ott.shapes.Sphere(1.0);
  index_relative = 1.2;
  particle = ott.scat.shapeforce.Shape(shape, index_relative);

  beam.position = [10;0;0];
  force = particle.force(beam);
  testCase.verifyEqual(force, [0;0;0], 'AbsTol', 1e-15, 'far away');
end

function testForcePositionArray(testCase)

  beam = ott.beam.abstract.Gaussian(1.0);
  
  shape = ott.shapes.Sphere(1.0);
  index_relative = 1.2;
  particle = ott.scat.shapeforce.Shape(shape, index_relative);

  force = particle.force(beam, 'position', randn(3, 5));
  testCase.verifyEqual(size(force), [3, 5]);

end
