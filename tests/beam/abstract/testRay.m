function tests = testRay
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  beam = ott.beam.abstract.Ray();
  
  testCase.verifyEqual(beam.origin, beam.position);
  testCase.verifyEqual(beam.direction, beam.rotation * [0;0;1]);
  
  h = figure();
  beam.visualise();
  close(h);

end

function testConvertBeam(testCase)

  abs_beam = ott.beam.abstract.Ray();

  % Beam casts

  beam = ott.beam.Beam(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Ray');
  verifyProperties(testCase, ?ott.beam.abstract.Ray, beam, abs_beam);

end