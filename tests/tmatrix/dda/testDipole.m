function tests = testDipole
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruction(testCase)

  xyz = [0;0;0];
  polarization = [0;0;1];
  beam = ott.tmatrix.dda.Dipole(xyz, polarization);

  testCase.verifyEqual(beam.locations, xyz, 'xyz');
  testCase.verifyEqual(beam.polarization, polarization, 'pol');
  testCase.verifyEqual(beam.xySymmetry, false, 'xy-sym');
  testCase.verifyEqual(beam.zRotSymmetry, 1, 'z-sym');
  testCase.verifyEqual(beam.rorder, 0, 'rorder');
  testCase.verifyEqual(beam.parity, 'even', 'parity');
  testCase.verifyEqual(beam.ndipoles, 1, 'ndipoles');
  testCase.verifyEqual(beam.nbeams, 1, 'nbeams');

end

function testFieldCoverage(testCase)

  xyz = [0;0;0];
  polarization = [0;0;1];
  beam = ott.tmatrix.dda.Dipole(xyz, polarization);

  rtp = randn(2, 5);

  E = beam.efarfield(rtp);
  testCase.verifySize(E.vrtp, [3, 5], 'sz efarfield');

  E = beam.hfarfield(rtp);
  testCase.verifySize(E.vrtp, [3, 5], 'sz hfarfield');

  dxyz = randn(3, 2);
  E = beam.efield(dxyz);
  testCase.verifySize(E.vrtp, [3, 2], 'sz efield');

  E = beam.hfield(dxyz);
  testCase.verifySize(E.vrtp, [3, 2], 'sz hfield');

end

function testLowMemory(testCase)

  xyz = [0;0;0];
  polarization = [0;0;1];
  beam = ott.tmatrix.dda.Dipole(xyz, polarization, ...
      'rorder', 2, 'zRotSymmetry', 4, 'xySymmetry', true, 'parity', 'odd');

  rtp = randn(2, 5);

  E1 = beam.efarfield(rtp);
  E2 = beam.efarfield(rtp, 'low_memory', true);

  testCase.verifyEqual(E2, E1, 'AbsTol', 1e-15, 'efarfield');

end


