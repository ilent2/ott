function tests = testBsc
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testPmNearfield(testCase)

  target = ott.bsc.Bsc([1;0;0;0;0;0;0;0], [0;0;0;0;0;0;0;0]);
  xyz = randn(3, 100);
  E = target.efield(xyz);

  pmbeam = ott.bsc.Pointmatch.FromNearfield(xyz, E.vxyz, 1:8);

  testCase.verifyEqual(pmbeam.a, target.a, 'AbsTol', 1e-15);
  testCase.verifyEqual(pmbeam.b, target.b, 'AbsTol', 1e-15);

end

function testPmFarfield(testCase)

  target = ott.bsc.Bsc([1;0;0;0;0;0;0;0], [1i;0;0;0;0;0;0;0]);
  rtp = rand(2, 100)*pi;
  E = target.efarfield(rtp);

  pmbeam = ott.bsc.Pointmatch.FromFarfield(rtp, E.vrtp, 1:8);

  testCase.verifyEqual(pmbeam.a, target.a, 'AbsTol', 1e-15);
  testCase.verifyEqual(pmbeam.b, target.b, 'AbsTol', 1e-15);

end

