function tests = testBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  angle = pi/4;
  beam = ott.beam.abstract.Bessel(angle);
  
  testCase.verifyEqual(beam.angle, angle);
  
  h = figure();
  ott.beam.Ray(beam).visualise()
  axis equal;
  close(h);
end

function testCasts(testCase)

  abs_beam = ott.beam.abstract.Bessel(pi/4);

  % Beam casts

  beam = ott.beam.Ray(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Ray');
  verifyProperties(testCase, ?ott.beam.abstract.Bessel, beam, abs_beam);

  % VSWF
  
  beam = ott.beam.vswf.Bsc(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.Bessel');
  verifyProperties(testCase, ?ott.beam.abstract.Bessel, beam, abs_beam);
  
  beam = ott.beam.vswf.Bessil(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.Bessel');
  verifyProperties(testCase, ?ott.beam.abstract.Bessel, beam, abs_beam);
  
  beam = ott.beam.vswf.BessilBasis(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.BesselBasis');
  verifyProperties(testCase, ?ott.beam.abstract.Bessel, beam, abs_beam);
end
