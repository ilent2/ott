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
  testCase.verifyEqual(dipole.location, dipole.position);
  
  h = figure();
  dipole.visualise('axis', 'y', 'field', 'Re(Ez)');
  close(h);

end

function testCasts(testCase)

  abs_beam = ott.beam.abstract.Dipole('polarization', [0;0;1]);

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
    
  % VSWF

  beam = ott.beam.vswf.Bsc(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.Dipole');
  verifyProperties(testCase, ?ott.beam.abstract.Dipole, ...
      beam, abs_beam);

  beam = ott.beam.vswf.Dipole(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.Dipole');
  verifyProperties(testCase, ?ott.beam.abstract.Dipole, ...
      beam, abs_beam);
end

