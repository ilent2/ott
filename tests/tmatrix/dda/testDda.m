function tests = testDda
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testFromShape(testCase)

  shape = ott.shape.Sphere(0.1);
  dda = ott.tmatrix.dda.Dda.FromShape(shape, 'index_relative', 1.2);

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
  alpha = diag([3, 2, 1]);
  dda = ott.tmatrix.dda.Dda(xyz, alpha);

  Einc = [1;0;0];
  dipoles = dda.solve(Einc);

  testCase.verifyInstanceOf(dipoles, 'ott.tmatrix.dda.Dipole', 'instance');
  testCase.verifyEqual(dipoles.polarization, [3; 0; 0], 'polarization');

end

function testTwoDipolesMirror(testCase)

  xyz1 = [0;0;1];
  xyz2 = [xyz1, -xyz1];
  alpha1 = eye(3);
  alpha2 = [alpha1, alpha1];
  
  dda1 = ott.tmatrix.dda.Dda(xyz1, alpha1, 'xySymmetry', true);
  dda2 = ott.tmatrix.dda.Dda(xyz2, alpha2);
  
  Einc = [1;0;1];
  
  dipoles1 = dda1.solve(Einc, 'parity', 'odd');
  dipoles2 = dda2.solve([Einc; [-1;-1;1].*Einc]);
  
  testCase.verifyEqual(dipoles1.polarization(1:3), ...
      dipoles2.polarization(1:3), 'AbsTol', 1e-16, 'e1');
  
  dipoles1 = dda1.solve(Einc, 'parity', 'even');
  dipoles2 = dda2.solve([Einc; [1;1;-1].*Einc]);
  
  testCase.verifyEqual(dipoles1.polarization(1:3), ...
      dipoles2.polarization(1:3), 'AbsTol', 1e-16, 'e2');

end

