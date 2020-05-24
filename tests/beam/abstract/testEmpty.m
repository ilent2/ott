function tests = testBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)
  beam = ott.beam.abstract.Empty();
end

function testCasts(testCase)

  beam = ott.beam.abstract.Empty();

  % Beam casts

  beam = ott.beam.Beam(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Empty');
  verifyProperties(testCase, ?ott.beam.abstract.Empty, ...
      beam, abs_beam);

  beam = ott.beam.Empty(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Empty');
  verifyProperties(testCase, ?ott.beam.abstract.Empty, ...
      beam, abs_beam);

  % Empty ArrayTypes

  beam = ott.beam.PlaneWave(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.PlaneWave');
  testCase.verifyEmpty(beam);
  verifyProperties(testCase, ?ott.beam.abstract.Empty, ...
      beam, abs_beam);

  beam = ott.beam.Bsc(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Bsc');
  testCase.verifyEmpty(beam);
  verifyProperties(testCase, ?ott.beam.abstract.Empty, ...
      beam, abs_beam);

  beam = ott.beam.Dipole(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Dipole');
  testCase.verifyEmpty(beam);
  verifyProperties(testCase, ?ott.beam.abstract.Empty, ...
      beam, abs_beam);

  beam = ott.beam.Ray(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Ray');
  testCase.verifyEmpty(beam);
  verifyProperties(testCase, ?ott.beam.abstract.Empty, ...
      beam, abs_beam);
end
