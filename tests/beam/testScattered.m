function tests = testScattered
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testFromParticle(testCase)

  shape = ott.shape.Sphere(1e-6);
  particle = ott.particle.Fixed.FromShape(shape, 'index_relative', 1.2);

  beam = ott.beam.Gaussian();
  sbeam = particle * beam;

  testCase.verifyInstanceOf(sbeam, 'ott.beam.Scattered');
end

