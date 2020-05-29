function tests = testPlaneWave
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  beam = ott.beam.abstract.PlaneWave();
  
  testCase.verifyEqual(beam.origin, beam.position);
  testCase.verifyEqual(beam.direction, beam.rotation * [0;0;1]);
  
  h = figure();
  beam.visualise('axis', 'y', 'field', 'Re(Ey)');
  close(h);

end

function testConvertBeam(testCase)

  abs_beam = ott.beam.abstract.PlaneWave();

  % Beam casts

  beam = ott.beam.Beam(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.PlaneWave');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);

end
