function tests = testArrays
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  beam1 = ott.beam.BscBeam();
  beam2 = ott.beam.BscBeam();
  beams = [beam1, beam2];
  
  beamA = ott.beam.Array(beams);
  
  testCase.verifyEqual(beamA.data, beams, 'data');
  testCase.verifyEqual(beamA.arrayType, 'array', 'type');
end

function testCoherent(testCase)

  beam(1) = ott.beam.Gaussian();
  beam(2) = ott.beam.Gaussian();
  beam(1).position = [-1;0;0]*1e-6;
  beam(2).position = [1;0;0]*1e-6;

  e = beam.visNearfield('range', [1,1]*5e-6, 'field', 'Ex');
  es = sum(e, 3);

  cbeam = ott.beam.Coherent(beam);
  e = cbeam.visNearfield('range', [1, 1]*5e-6, 'field', 'Ex');

  testCase.verifyEqual(e, es, 'Ex');

end

function testIncoherent(testCase)

  beam(1) = ott.beam.Gaussian();
  beam(2) = ott.beam.Gaussian();
  beam(1).position = [-1;0;0]*1e-6;
  beam(2).position = [1;0;0]*1e-6;

  e = beam.visNearfield('range', [1,1]*5e-6, 'field', 'irradiance');
  es = sum(e, 3);

  cbeam = ott.beam.Incoherent(beam);
  e = cbeam.visNearfield('range', [1, 1]*5e-6, 'field', 'irradiance');

  testCase.verifyEqual(e, es, 'irradiance');

end

function testPosition(testCase)

  offset = [0;1;0]*1e-6;

  beam(1) = ott.beam.Gaussian();
  beam(2) = ott.beam.Gaussian();
  beam(1).position = [-1;0;0]*1e-6;
  beam(2).position = [1;0;0]*1e-6;
  cbeam1 = ott.beam.Incoherent(beam);
  cbeam1.position = offset;

  beam = beam.translateXyz(offset);
  cbeam2 = ott.beam.Incoherent(beam);

  e1 = cbeam1.visNearfield();
  e2 = cbeam2.visNearfield();

  testCase.verifyEqual(e1, e2, 'position incorrect');

end

