function tests = testBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  beam = ott.beam.abstract.PlaneWave();
  
  testCase.verifyEqual(beam.origin, beam.position);
  testCase.verifyEqual(beam.direction, beam.rotation * [0;0;1]);

end

function testConvertBeam(testCase)

  abs_beam = ott.beam.abstract.PlaneWave();

  % Beam casts

  beam = ott.beam.Beam(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.PlaneWave');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);

  beam = ott.beam.PlaneWave(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.PlaneWave');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);

  beam = ott.beam.Ray(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Ray');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);

  % VSWF Casts

  beam = ott.beam.vswf.Bsc(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.PlaneWave');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);

  beam = ott.beam.vswf.PlaneWave(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.PlaneWave');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);

  beam = ott.beam.vswf.PlaneBasis(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.PlaneBasis');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);

end
