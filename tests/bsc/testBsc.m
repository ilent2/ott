function tests = testBsc
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstructEmpty(testCase)

  bsc = ott.bsc.Bsc();
  testCase.verifyEmpty(bsc.a, 'invalid a');
  testCase.verifyEmpty(bsc.b, 'invalid b');
  testCase.verifyEqual(bsc.Nmax, 0, 'invalid Nmax');
  testCase.verifyEqual(bsc.power, 0, 'invalid power');

end

function testConstructArray(testCase)

  a = rand(3, 5);
  b = rand(3, 5);
  bsc = ott.bsc.Bsc(a, b);
  testCase.verifySize(bsc, [1, 5], 'bsc size');
  testCase.verifyEqual([bsc.a], a, 'a');
  testCase.verifyEqual([bsc.b], b, 'b');

end

function testPmNearfield(testCase)

  target = ott.bsc.Bsc([1;0;0;0;0;0;0;0], [0;0;0;0;0;0;0;0]);
  rtp = randn(3, 100);
  E = target.efieldRtp(rtp);

  pmbeam = full(ott.bsc.Bsc.PmNearfield(rtp, E.vrtp, 1:8));

  testCase.verifyEqual(pmbeam.a, target.a, 'AbsTol', 1e-15);
  testCase.verifyEqual(pmbeam.b, target.b, 'AbsTol', 1e-15);

end

function testPmFarfield(testCase)

  target = ott.bsc.Bsc([1;0;0;0;0;0;0;0], [1i;0;0;0;0;0;0;0]);
  rtp = rand(2, 100)*pi;
  E = target.efarfield(rtp);

  pmbeam = full(ott.bsc.Bsc.PmFarfield(rtp, E.vrtp, 1:8));

  testCase.verifyEqual(pmbeam.a, target.a, 'AbsTol', 1e-15);
  testCase.verifyEqual(pmbeam.b, target.b, 'AbsTol', 1e-15);

end

function testConstructCopy(testCase)

  bsc = ott.bsc.Bsc([1;2;3], [0;0;0]);
  bscCopy = ott.bsc.Bsc(bsc);
  
  testCase.verifyEqual(bscCopy, bsc);
end

function testConstructFromDense(testCase)

  a = [1, 2, 3, 4, 5, 7, 8].';
  b = 10+a;
  n = [1, 1, 1, 2, 2, 2, 2].';
  m = [-1, 0, 1, -2, -1, 1, 2].';
  bsc = ott.bsc.Bsc.FromDenseBeamVectors(a, b, n, m);
  testCase.verifyEqual(bsc.a, sparse([1, 2, 3, 4, 5, 0, 7, 8].'), 'invalid a');
  testCase.verifyEqual(bsc.b, sparse(10+[1, 2, 3, 4, 5, -10, 7, 8].'), 'invalid b');
  testCase.verifyEqual(bsc.Nmax, 2, 'invalid Nmax');
  testCase.verifyEqual(bsc.power, sum(a.^2 + b.^2), 'invalid power');

  [on, om] = bsc.getModeIndices();
  testCase.verifyEqual(on, n, 'n');
  testCase.verifyEqual(om, m, 'm');

end

function testFromDenseEmpty(testCase)

  a = [];
  b = [];
  nn = [];
  mm = [];

  bsc = ott.bsc.Bsc.FromDenseBeamVectors(a, b, nn, mm);
  testCase.verifyEmpty(bsc.a);
  testCase.verifyEmpty(bsc.b);
  testCase.verifyEqual(bsc.Nmax, 0);

end

function testBasisSet(testCase)

  bsc = ott.bsc.Bsc.BasisSet(1:3);
  target = ott.bsc.Bsc([eye(3), zeros(3)], [zeros(3), eye(3)]);
  
  testCase.verifyEqual(bsc, sparse(target));
end

function testSetNmax(testCase)

  bsc = ott.bsc.Bsc([10;10;10;0.1]);

  bsc.Nmax = 20;
  testCase.verifyEqual(bsc.Nmax, 20);

  bsc.Nmax = 10;
  testCase.verifyEqual(bsc.Nmax, 10);

end

function testShrinkNmax(testCase)

  beam = ott.bsc.Bsc([10;10;10;0.1]);
  testCase.assertEqual(beam.Nmax, 2);
  
  sbeam = beam.shrinkNmax('RelTol', 1.0e-2);
  testCase.verifyEqual(sbeam.Nmax, 1);
  
  sbeam = beam.shrinkNmax('AbsTol', 1, 'RelTol', []);
  testCase.verifyEqual(sbeam.Nmax, 1);

end

function testConcateration(testCase)

  a = ott.bsc.Bsc([1, 2, 3].', [4, 5, 6].');
  b = ott.bsc.Bsc([1, 2, 3, 7, 8, 9, 10, 11].', [4, 5, 6].');
  bsc = [a, b];

  testCase.verifySize(bsc, [1, 2]);
  testCase.verifyEqual(bsc(1), a);
  testCase.verifyEqual(bsc(2), b);

end

function testForceTorque(testCase)
  % Verify the result from forcetorque matches individual force/torque funcs

  a = ott.bsc.Bsc([1; 2; 3], [4; 5; 6]);
  b = ott.bsc.Bsc([1; 2; 3; 7; 8; 9; 10; 11], [4; 5; 6]);

  f0 = a.force(b);
  t0 = a.torque(b);
  s0 = a.spin(b);

end

function testGetCoefficients(testCase)

  abeam = ott.bsc.Bsc([1; 2; 3], [4; 5; 6]);
  bbeam = ott.bsc.Bsc([1; 2; 3; 7; 8; 9; 10; 11], [4; 5; 6]);
  beam = [abeam, bbeam];

  ab = beam.getCoefficients([1, 2, 3]);
  testCase.verifyEqual(ab, [1, 1; 2, 2; 3, 3; 4, 4; 5, 5; 6, 6], 'ab');

  [a, b] = beam.getCoefficients([1, 2, 3]);
  testCase.verifyEqual(a, [1, 1; 2, 2; 3, 3], 'a');
  testCase.verifyEqual(b, [4, 4; 5, 5; 6, 6], 'b');
  
  [a, b] = beam.getCoefficients(3);
  testCase.verifyEqual(a, [3, 3], 'a');
  testCase.verifyEqual(b, [6, 6], 'b');

  ab = beam.getCoefficients();
  testCase.verifySize(ab, [16, 2], 'size');

end

function testPlus(testCase)

  abeam = ott.bsc.Bsc([1; 2; 3]);
  bbeam = ott.bsc.Bsc([4; 5; 6]);

  sbeam = abeam + bbeam;
  testCase.verifyEqual(sbeam.a, [5; 7; 9], 'a only');

  sbeam = [abeam, abeam] + bbeam;
  testCase.verifyEqual([sbeam.a], [5, 5; 7, 7; 9, 9], 'a only in array');

end

function testMinus(testCase)

  abeam = ott.bsc.Bsc([1; 2; 3]);
  bbeam = ott.bsc.Bsc([4; 5; 6]);

  sbeam = abeam - bbeam;
  testCase.verifyEqual(sbeam.a, [-3; -3; -3], 'a only');

end

function testDivide(testCase)

  abeam = ott.bsc.Bsc([2; 4; 8]);

  rbeam = abeam ./ 2;
  testCase.verifyEqual(rbeam.a, [1;2;4]);

  rbeam = abeam / 2;
  testCase.verifyEqual(rbeam.a, [1;2;4]);

end

function testGpuArray(testCase)

  abeam = ott.bsc.Bsc([1; 2; 3]);
  abeam = gpuArray(abeam);
  testCase.verifyInstanceOf(abeam.a, 'gpuArray', 'cast to gpu');
  abeam = gather(abeam);
  testCase.verifyInstanceOf(abeam.a, 'double', 'cast from gpu');

end

function testSparse(testCase)

  beam = ott.bsc.Bsc([1; 2; 3]);
  testCase.verifyFalse(issparse(beam), 'full check');

  sbeam = sparse(beam);
  testCase.verifyTrue(issparse(sbeam), 'sparse check');
  testCase.verifyTrue(issparse(sbeam.a), 'is sparse');
  testCase.verifyEqual(full(sbeam.a), beam.a, 'sparse');

  fbeam = full(sbeam);
  testCase.verifyFalse(issparse(fbeam), 'after cast check');

  sbeam2 = beam.makeSparse();
  testCase.verifyTrue(issparse(sbeam2), 'sparse check');

  % Again with beam array

  beam = [beam, beam];
  testCase.verifyFalse(any(issparse(beam)), 'A: full check');

  beam2 = sparse(beam);
  testCase.verifyTrue(all(issparse(beam2)), 'A: sparse check');

  beam3 = full(beam2);
  testCase.verifyFalse(any(issparse(beam3)), 'A: full check after');

  beam4 = beam.makeSparse();
  testCase.verifyTrue(all(issparse(beam4)), 'A: make sparse');

end

function testSum(testCase)

  Nmax = 5;
  ab1 = randn(ott.utils.combined_index(Nmax, Nmax), 1);
  ab2 = randn(ott.utils.combined_index(Nmax, Nmax), 1);
  beam1 = ott.bsc.Bsc(ab1, ab1);
  beam2 = ott.bsc.Bsc(ab2, ab2);

  beam3 = beam1 + beam2;
  testCase.assertSize(beam3, [1, 1], 'sz');

  beamU = [beam1, beam2];
  testCase.verifyEqual(beamU.sum(), beam3, ...
    'beamU.sum() incorrect');
  testCase.verifyEqual(sum(beamU), beam3, ...
    'sum(beamU) incorrect');

end

function testPower(testCase)

  abeam = ott.bsc.Bsc([1; 2; 3]);
  abeam.power = 1;
  testCase.verifyEqual(abeam.power, 1);

end

function testTranslateXyz(testCase)

  abeam = ott.bsc.Bsc([1; 2; 3]);
  bbeam = abeam.translateXyz([0;0;1]);

end

function testRealImagAbs(testCase)

  beam = ott.bsc.Bsc([1i; 2; 3i]);

  abeam = abs(beam);
  testCase.verifyEqual(abeam.a, [1; 2; 3], 'abs');

  abeam = real(beam);
  testCase.verifyEqual(abeam.a, [0; 2; 0], 'real');

  abeam = imag(beam);
  testCase.verifyEqual(abeam.a, [1; 0; 3], 'imag');

end

function testTmatrixCast(testCase)

  beam = ott.bsc.Bsc([1i; 2; 3i]);
  beam = [beam, beam, beam, beam, beam, beam];

  tmatrix = ott.tmatrix.Tmatrix(beam);
  testCase.verifyEqual(tmatrix.Nmax, [1, 1], 'Nmax');

end

function testHfieldRtp(testCase)

  beam = ott.bsc.Bsc([1i; 2; 3i]);
  rtp = randn(3, 5);
  E = beam.hfieldRtp(rtp);
  
end

function testDipoleFieldValues(testCase)

  beams = ott.bsc.Bsc.BasisSet(1:3);
  E0 = beams.efieldRtp([0;0;0]);
  
  I = vecnorm(E0.vxyz(1:3, :));
  
  % Compare to dipole field values from OTTv1.5
  testCase.verifyEqual(I, [0, 0, 0, 0.2303, 0.2303, 0.2303], ...
      'AbsTol', 1e-4, 'dipole field values');

end
