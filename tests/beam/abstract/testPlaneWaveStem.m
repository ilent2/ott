function tests = testPlaneWaveStem
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConvertBeam(testCase)

  abs_beam = ott.beam.abstract.PlaneWave();

  % Beam casts

  beam = ott.beam.Ray(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Ray');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);

  beam = ott.beam.PlaneWave(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.PlaneWave');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);

  % Array types (should match ott.beam.Beam cast)
  
  beam = ott.beam.Array(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.PlaneWave');
  testCase.verifyEqual(beam.array_type, 'array', 'atype');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);

  beam = ott.beam.Incoherent(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.PlaneWave');
  testCase.verifyEqual(beam.array_type, 'incoherent', 'atype');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);

  beam = ott.beam.Coherent(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.PlaneWave');
  testCase.verifyEqual(beam.array_type, 'coherent', 'atype');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);
  
  % VSWF
  
  beam = ott.beam.vswf.Bsc(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.PlaneWave');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);
  
  beam = ott.beam.vswf.PlaneWave(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.PlaneWave');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);
  
  beam = ott.beam.vswf.PlaneBasis(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.PlaneBasis');
  
  % Abstract
  
  beam = ott.beam.vswf.Ray(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.abstract.Ray');
  verifyProperties(testCase, ?ott.beam.abstract.PlaneWave, beam, abs_beam);
  
  abs_beam = ott.beam.abstract.Ray();
  beam = ott.beam.vswf.Ray(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.abstract.PlaneWave');
  verifyProperties(testCase, ?ott.beam.abstract.Ray, beam, abs_beam);

end