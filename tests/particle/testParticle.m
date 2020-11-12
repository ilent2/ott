function tests = testParticle
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');

  shape = ott.shape.Sphere(1e-6);
  index_relative = 1.2;
  mass = 1.7;
  testCase.TestData.fixed = ott.particle.Fixed.FromShape(shape, ...
      'index_relative', index_relative, 'internal', true, ...
      'mass', mass, 'wavelength0', 1064e-9);
    
end

function testSurf(testCase)

  f = figure();
  testCase.addTeardown(@() close(f));
  
  testCase.TestData.fixed.surf();

end

function testSetMassFromDensity(testCase)

  particle = testCase.TestData.fixed;
  
  particle.shape.volume = 1;
  density = 5;
  particle = particle.setMassFromDensity(density);
  testCase.verifyEqual(particle.mass, 5, 'AbsTol', 1e-14, 'mass');

end
