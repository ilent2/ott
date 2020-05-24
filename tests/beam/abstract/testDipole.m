function tests = testBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  polarization = [0;0;1];
  dipole = ott.beam.abstract.Dipole('polarization', polarization);
  testCase.verifyEqual(dipole.polarization, polarization);

end

function testCasts(testCase)

  dipole = ott.beam.abstract.Dipole('polarization', [0;0;1]);

  beam = ott.beam.Beam(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Dipole');
  verifyProperties(testCase, ?ott.beam.abstract.Dipole, ...
      beam, abs_beam);

  beam = ott.beam.Dipole(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Dipole');
  verifyProperties(testCase, ?ott.beam.abstract.Dipole, ...
      beam, abs_beam);

  beam = ott.beam.Coherent(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.Dipole');
  verifyProperties(testCase, ?ott.beam.abstract.Dipole, ...
      beam, abs_beam);

  beam = ott.beam.vswf.Bsc(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.Bsc');
  verifyProperties(testCase, ?ott.beam.abstract.Dipole, ...
      beam, abs_beam);
end

