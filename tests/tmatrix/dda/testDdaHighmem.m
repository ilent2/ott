function tests = testDdaHighmem
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testFromShape(testCase)

  shape = ott.shape.Sphere(0.1);
  dda = ott.tmatrix.dda.DdaHighMem.FromShape(shape, 'relative_index', 1.2);

  testCase.verifyEqual(dda.zRotSymmetry, 0, 'zrot');
  testCase.verifyEqual(dda.xySymmetry, true, 'xysym');

  ndp = dda.ndipoles;
  testCase.verifySize(dda.interaction, [3*ndp, 3*ndp, 4, 2], 'sz int');

end

function testSolve(testCase)

  xyz = [0;0;0];
  alpha = eye(3);
  dda1 = ott.tmatrix.dda.Dda(xyz, alpha);
  dda2 = ott.tmatrix.dda.DdaHighMem(xyz, alpha);

  beam = ott.bsc.Bsc(1, 0);
  Einc = beam.efield(xyz);
  dipoles1 = dda1.solve(Einc);
  dipoles2 = dda2.solve(Einc);

  testCase.verifyEqual(dipoles1, dipoles2, 'match low memory');

end

