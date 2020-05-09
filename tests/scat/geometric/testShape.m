function tests = testShape
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruct(testCase)

  shape = ott.shapes.Sphere(1.0);
  index_relative = 1.5;
  
  particle = ott.scat.geometric.Shape(shape, index_relative);
  testCase.verifyEqual(particle.shape, shape, 'shape');
  testCase.verifyEqual(particle.index_relative, index_relative);
end

function testScatter(testCase)

  shape = ott.shapes.Sphere(1.0);
  index_relative = 1.5;
  particle = ott.scat.geometric.Shape(shape, index_relative);
  beam = ott.beam.Ray('direction', [0;0;1], 'origin', [0;0;-3]);
  
  [sbeam, ibeam] = particle.scatter(beam, 'max_iterations', 20);
  
  testCase.verifyEqual(sbeam.power, beam.power, ...
    'AbsTol', 2.0e-2, 'power prserved');
  testCase.verifyEqual(size(ibeam), [0, 0], 'no internal rays');
  
end

function testForceTorque(testCase)

  shape = ott.shapes.Sphere(1.0);
  index_relative = 1.5;
  particle = ott.scat.geometric.Shape(shape, index_relative);
  beam = ott.beam.Ray('direction', [0;0;1], 'origin', [0;0;-3]);

  force1 = particle.force(beam);
  testCase.verifyEqual(size(force1), [3, 1], 'size');
  testCase.verifyEqual(force1(1:2), [0;0], ...
    'AbsTol', 1e-15, 'off-axis terms');
  
  force2 = beam.force(particle);
  testCase.verifyEqual(force2, -force1, ...
    'AbsTol', 1e-15, 'oposite sign');
  
  torque = particle.torque(beam);
  testCase.verifyEqual(torque, zeros(3, 1), ...
    'AbsTol', 1e-15, 'torque');
  
  [force3, torque3] = particle.forcetorque(beam);
  testCase.verifyEqual(force3, force1, ...
    'AbsTol', 1e-15, 'ft force');
  testCase.verifyEqual(torque3, zeros(3, 1), ...
    'AbsTol', 1e-15, 'ft torque');
end

function testForceMie(testCase)

  % TODO: We might want to choose a better beam for comparison

  shape = ott.shapes.Sphere(1.0);
  index_relative = 1.5;
  particleGo = ott.scat.geometric.Shape(shape, index_relative);
  particleMie = ott.scat.vswf.Mie(shape, 'index_relative', index_relative);
  beam = ott.beam.PlaneWave('direction', [0;0;1], 'origin', [0;0;-3]);
  
  force1 = particleGo.force(beam);
  force2 = particleMie.force(beam);
  testCase.verifyEqual(force1, force2, 'force');
  
  torque1 = particleGo.torque(beam);
  torque2 = particleMie.torque(beam);
  testCase.verifyEqual(torque1, torque2, 'torque');
end
