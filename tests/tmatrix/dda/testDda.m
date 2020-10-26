function tests = testDda
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testFromShape(testCase)

  shape = ott.shape.Sphere(0.1);
  dda = ott.tmatrix.dda.Dda.FromShape(shape, 'relative_index', 1.2);

  testCase.verifyEqual(dda.zRotSymmetry, 0, 'zrot');
  testCase.verifyEqual(dda.xySymmetry, true, 'xysym');

end

function testConstruct(testCase)

  xyz = [1;0;0];
  alpha = randn(3, 3);
  dda = ott.tmatrix.dda.Dda(xyz, alpha, 'zRotSymmetry', 2);

  testCase.verifyEqual(dda.ndipoles, 1, 'ndipoles');
  testCase.verifyEqual(dda.locations, xyz, 'xyz');
  testCase.verifyEqual(dda.polarizability, alpha, 'alpha');
  testCase.verifyEqual(dda.zRotSymmetry, 2, 'zrot');

end

function testSolve(testCase)

  xyz = [0;0;0];
  alpha = eye(3);
  dda = ott.tmatrix.dda.Dda(xyz, alpha);

  beam = ott.bsc.Bsc(1, 0);
  Einc = beam.efield(xyz);
  dipoles = dda.solve(Einc);

  testCase.verifyInstance(dipoles, 'ott.tmatrix.dda.Dipole', 'instance');

end

