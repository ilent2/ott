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
  beam = beam.shrinkNmax('RelTol', 1.0e-2);
  testCase.verifyEqual(beam.Nmax, 1);

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

  [f, t, s] = a.forcetorque(b);
  f0 = a.force(b);
  t0 = a.torque(b);
  s0 = a.spin(b);
  testCase.verifyEqual(f0, f, 'force');
  testCase.verifyEqual(t0, t, 'torque');
  testCase.verifyEqual(s0, s, 'spin');

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

  abeam = ott.bsc.Bsc([1; 2; 3]);
  testCase.verifyFalse(issparse(abeam), 'full check');

  sbeam = sparse(abeam);
  testCase.verifyTrue(issparse(sbeam), 'sparse check');
  testCase.verifyTrue(issparse(sbeam.a), 'is sparse');
  testCase.verifyEqual(full(sbeam.a), abeam.a, 'sparse');

  fbeam = full(sbeam);
  testCase.verifyFalse(issparse(fbeam), 'after cast check');

  sbeam2 = abeam.makeSparse();
  testCase.verifyTrue(issparse(sbeam2), 'sparse check');

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

function testVisNearfield(testCase)

  f = figure();
  testCase.addTeardown(@close, f);

  abeam = ott.bsc.Bsc([1; 2; 3]);
  abeam.visNearfield();

  abeam.visNearfield('axis', 'x', 'range', [5, 5]);
  abeam.visNearfield('axis', {[1;0;0], [0;1;0]}, 'range', [-1, 1, -1, 1]);
  abeam.visNearfield('axis', {[1;0;0], [0;1;0], [0;0;1]}, ...
      'range', {linspace(0, 1), linspace(0, 2)});

end

function testVisFarfield(testCase)

  f = figure();
  testCase.addTeardown(@close, f);

  abeam = ott.bsc.Bsc([1; 2; 3]);
  abeam.visFarfield();

end

function testVisFarfieldSlice(testCase)

  f = figure();
  testCase.addTeardown(@close, f);

  abeam = ott.bsc.Bsc([1; 2; 3]);
  abeam.visFarfieldSlice(0.0);

end

function testVisFarfieldSphere(testCase)

  f = figure();
  testCase.addTeardown(@close, f);

  abeam = ott.bsc.Bsc([1; 2; 3]);
  abeam.visFarfieldSphere();

end

function testPower(testCase)

  abeam = ott.bsc.Bsc([1; 2; 3]);
  abeam.power = 1;
  testCase.verifyEqual(abeam.power, 1);

end

function testIntensityMoment(testCase)

  abeam = ott.bsc.Bsc([1; 2; 3]);
  [moment, ints] = abeam.intensityMoment();

end

function testTranslateXyz(testCase)

  abeam = ott.bsc.Bsc([1; 2; 3]);
  bbeam = abeam.translateXyz([0;0;1]);

end

function testFieldCoverage(testCase)

  abeam = ott.bsc.Bsc([1; 2; 3]);
  H = abeam.hfield([0;0;1]);
  [E, H] = abeam.ehfield([0;0;1]);
  [Ep, Hp] = abeam.ehparaxial([0;0.5]);
  Hp = abeam.hparaxial([0;0.5]);
  [Ef, Hf] = abeam.ehfarfield([0;0.5]);
  Hf = abeam.hfarfield([0;0.5]);

end

