function tests = testBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  radius = 1.0;
  beam = ott.beam.abstract.TopHat(radius);
  
  testCase.verifyEqual(beam.radius, radius);
  testCase.verifyEqual(beam.power, pi, 'power');
  
  h = figure();
  ott.beam.Ray(beam).visualise()
  axis equal;
  close(h);

end

function testCasts(testCase)

  abs_beam = ott.beam.abstract.TopHat(1.0);

  % Beam casts

  beam = ott.beam.Beam(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.abstract.MaskedNearfield3d');
  verifyProperties(testCase, ?ott.beam.abstract.Beam, beam, abs_beam);

  beam = ott.beam.Ray(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Ray');
  verifyProperties(testCase, ?ott.beam.abstract.Beam, beam, abs_beam);

  % VSWF
  
  beam = ott.beam.vswf.Pointmatch(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.NearfieldPm');
  verifyProperties(testCase, ?ott.beam.abstract.Beam, beam, abs_beam);
  
  beam = ott.beam.vswf.NearfieldPm(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.NearfieldPm');
  verifyProperties(testCase, ?ott.beam.abstract.Beam, beam, abs_beam);
  
  % Abstract
  
  beam = ott.beam.abstract.MaskedNearfield3d(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.abstract.MaskedNearfield3d');
  verifyProperties(testCase, ?ott.beam.abstract.Beam, beam, abs_beam);

end

