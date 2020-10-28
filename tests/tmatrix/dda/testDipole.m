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

function testFieldSize(testCase)

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
  polarization = [1;2;3];
  beam = ott.tmatrix.dda.Dipole(xyz, polarization, ...
      'rorder', 2, 'zRotSymmetry', 4, 'xySymmetry', true, 'parity', 'odd');

  rtp = randn(2, 5);

  E1 = beam.efarfield(rtp, 'low_memory', false);
  E2 = beam.efarfield(rtp, 'low_memory', true);

  testCase.verifyEqual(E2, E1, 'AbsTol', 1e-15, 'efarfield');

end

function testDipoleFarfield(testCase)

  distance = 1e6;
  target_rtp = [distance; 0; 0];
  target_xyz = ott.utils.rtp2xyz(target_rtp);

  xyz = [0;0;1];
  polarization = [1;1;1];
  beam = ott.tmatrix.dda.Dipole(xyz, polarization);

  % Phase factor for near-field
  P = exp(1i*2*pi*distance)./distance;

  farfield = beam.efarfield(target_rtp);
  nearfield = beam.efield(target_xyz);
  testCase.verifyEqual(farfield.vxyz, nearfield./P, 'AbsTol', 1.0e-3, 'e');

  farfield = beam.hfarfield(target_rtp);
  nearfield = beam.hfield(target_xyz);
  testCase.verifyEqual(farfield.vxyz, nearfield./P, 'AbsTol', 1.0e-3, 'h');

end
