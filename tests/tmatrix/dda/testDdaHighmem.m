function tests = testDdaHighmem
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testFromShape(testCase)

  shape = ott.shape.Sphere(0.1);
  dda = ott.tmatrix.dda.DdaHighMem.FromShape(shape, 'index_relative', 1.2);

  testCase.verifyEqual(dda.zRotSymmetry, 0, 'zrot');
  testCase.verifyEqual(dda.xySymmetry, true, 'xysym');

  ndp = dda.ndipoles;
  testCase.verifySize(dda.interaction, [3*ndp, 3*ndp, 4, 2], 'sz int');

end

function testSolve(testCase)

  xyz = [1;1;1];
  alpha = eye(3);
  dda1 = ott.tmatrix.dda.Dda(xyz, alpha, ...
      'xySymmetry', true, 'zRotSymmetry', 4);
  dda2 = ott.tmatrix.dda.DdaHighMem(xyz, alpha, ...
      'xySymmetry', true, 'zRotSymmetry', 4);

  Einc = [1;2;3];
  dipoles1 = dda1.solve(Einc);
  dipoles2 = dda2.solve(Einc);

  testCase.verifyEqual(dipoles1.polarization, ...
      dipoles2.polarization, 'AbsTol', 1e-16, 'match low memory');

end

