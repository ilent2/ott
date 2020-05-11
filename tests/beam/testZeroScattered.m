function tests = testZeroScattered
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstructor(testCase)

  sphere = ott.shapes.Sphere(1.0);
  index_relative = 1.33;
  type = 'total';
  ibeam = ott.beam.PlaneWave();
  particle = ott.scat.shapeforce.Shape(sphere, index_relative);
  
  beam = ott.beam.ZeroScattered(type, ibeam, particle);

  testCase.verifyEqual(beam.particle, particle, 'particle');
  testCase.verifyEqual(beam.power, ibeam.power, 'power');
end

function testForce(testCase)

  sphere = ott.shapes.Sphere(1.0);
  index_relative = 1.33;
  type = 'total';
  ibeam = ott.beam.PlaneWave();
  particle = ott.scat.shapeforce.Shape(sphere, index_relative);
  beam = ott.beam.ZeroScattered(type, ibeam, particle);
  
  F = beam.force();
  targetF = particle.force(ibeam);
  testCase.verifyEqual(F, -targetF, 'force');
  
  F = beam.torque();
  targetF = particle.torque(ibeam);
  testCase.verifyEqual(F, -targetF, 'torque');
  
  F = beam.forcetorque();
  targetF = particle.forcetorque(ibeam);
  testCase.verifyEqual(F, -targetF, 'forcetorque');
end

function testVisualise(testCase)

  sphere = ott.shapes.Sphere(1.0);
  index_relative = 1.33;
  type = 'total';
  ibeam = ott.beam.PlaneWave();
  particle = ott.scat.shapeforce.Shape(sphere, index_relative);
  beam = ott.beam.ZeroScattered(type, ibeam, particle);
  
  im = beam.visualise();
  imTarget = ibeam.visualise();
  testCase.verifyEqual(im, imTarget, 'visualise');
end

