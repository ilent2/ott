function tests = testBsc
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testSmartCylinder(testCase)

  shape = ott.shape.Cylinder(0.5, 0.5);
  tmatrix = ott.tmatrix.Tmatrix.SmartCylinder(shape, 1.2);
  
  tmatrix = ott.tmatrix.Tmatrix.SmartCylinder(shape, 1.2, ...
      'tolerance', 'one');

end

function testFromShape(testCase)

  shape = ott.shape.Sphere(1.0);
  tmatrix = ott.tmatrix.Tmatrix.FromShape(shape, 1.2);
  testCase.verifyInstanceOf(tmatrix, 'ott.tmatrix.Mie', 'mie');
  
  shape = ott.shape.Cube(0.1);
  tmatrix = ott.tmatrix.Tmatrix.FromShape(shape, 1.2);
  testCase.verifyInstanceOf(tmatrix, 'ott.tmatrix.Pointmatch', 'cube');
end

function testConstructEmpty(testCase)

  tmatrix = ott.tmatrix.Tmatrix();
  testCase.verifyEmpty(tmatrix.data, 'invalid a');
  testCase.verifyEqual(tmatrix.Nmax, [0, 0], 'invalid Nmax');
  testCase.verifyEqual(tmatrix.type, 'scattered', 'invalid type');

end

function testConstructArray(testCase)

  data = randn(6);
  tmatrix = ott.tmatrix.Tmatrix({data, []});
  testCase.verifyEqual(tmatrix(1).Nmax, [1,1], 'nmax 1');
  testCase.verifyEqual(tmatrix(2).Nmax, [0,0], 'nmax 2');

end

function testBscCast(testCase)

  data = rand(16, 6);
  tmatrix = ott.tmatrix.Tmatrix(data);
  beam = ott.bsc.Bsc(tmatrix);

  testCase.verifySize(beam, [1, 6], 'size');
  testCase.verifyEqual(beam.getCoefficients(), data, 'data');

end

function testDowncast(testCase)

  data = rand(16, 6);
  tmatrix = ott.tmatrix.Tmatrix(data, 'type', 'total');
  tmatrix2 = ott.tmatrix.Tmatrix(tmatrix);

  testCase.verifyEqual(tmatrix2, tmatrix);

end

function testSparse(testCase)

  data = rand(16, 6);
  tmatrix = ott.tmatrix.Tmatrix(data, 'type', 'total');
  testCase.verifyFalse(issparse(tmatrix), 'full check');

  tmatrix2 = sparse(tmatrix);
  testCase.verifyTrue(issparse(tmatrix2), 'sparse check');
  testCase.verifyTrue(issparse(tmatrix2.data), 'is sparse');
  testCase.verifyEqual(full(tmatrix2.data), tmatrix.data, 'sparse');

  tmatrix3 = full(tmatrix2);
  testCase.verifyFalse(issparse(tmatrix3), 'after cast check');

  tmatrix4 = tmatrix.makeSparse();
  testCase.verifyTrue(issparse(tmatrix4), 'sparse check');

end

function testGpuArray(testCase)

  data = rand(16, 6);
  tmatrix = ott.tmatrix.Tmatrix(data, 'type', 'total');
  tmatrix = gpuArray(tmatrix);
  testCase.verifyInstanceOf(tmatrix.data, 'gpuArray', 'cast to gpu');
  tmatrix = gather(tmatrix);
  testCase.verifyInstanceOf(tmatrix.data, 'double', 'cast from gpu');

end

function testScatteredToTotal(testCase)

  tmatrix = ott.tmatrix.Tmatrix([], 'type', 'scattered');
  testCase.assertEqual(tmatrix.Nmax, [0,0], 'original Nmax');
  testCase.assertEqual(tmatrix.type, 'scattered');

  total = tmatrix.total;
  testCase.verifyEqual(total.data, [], 'empty');
  testCase.verifyEqual(total.type, 'total', 'new type');

  total.Nmax = 1;
  testCase.verifyEqual(total.data, eye(6), 'Nmax = 1');

end

function testTotalToScattered(testCase)

  tmatrix = ott.tmatrix.Tmatrix(eye(6), 'type', 'total');
  testCase.assertEqual(tmatrix.Nmax, [1,1], 'original Nmax');
  testCase.assertEqual(tmatrix.type, 'total');

  scattered = tmatrix.scattered;
  testCase.verifyEqual(scattered.Nmax, [1,1], 'Nmax');
  testCase.verifyEqual(scattered.data, zeros(6), 'data');
  testCase.verifyEqual(scattered.type, 'scattered', 'new type');

  scattered.Nmax = 0;
  testCase.verifyEqual(scattered.data, [], 'empty');

end

function testRealImagAbsFunctions(testCase)

  data = randn(6, 6) + 1i*randn(6, 6);
  T = ott.tmatrix.Tmatrix(data);

  Treal = real(T);
  testCase.verifyEqual(Treal.data, real(data), 'real');

  Timag = imag(T);
  testCase.verifyEqual(Timag.data, imag(data), 'imag');

  Tabs = abs(T);
  testCase.verifyEqual(Tabs.data, abs(data), 'abs');

end

function testShrinkPowerWarning(testCase)
  % Check shrinking gives a power warning

  T = ott.tmatrix.Tmatrix(eye(16), 'type', 'scattered');
  
  % Suppress no return warning
  warnState = warning('off', 'ott:utils:nargoutCheck:no_outputs');
  testCase.addTeardown(@warning, warnState);

  testCase.verifyWarning(@() T.setNmax(1), ...
      'ott:tmatrix:Tmatrix:setNmaxWarning:reltol', 'warn');

  testCase.verifyError(@() T.setNmax(1, 'powerloss', 'error'), ...
      'ott:tmatrix:Tmatrix:setNmaxWarning:reltol', 'error');

  testCase.verifyWarningFree(@() T.setNmax(1, 'powerloss', 'ignore'), ...
      'ignore');
end

function testPowerOnResize(testCase)

  % Check columnCheck
  T = ott.tmatrix.Tmatrix(eye(6), 'type', 'total');
  testCase.assertEqual(T.columnCheck(), ones(1, 6), 'AbsTol', 1.0e-15);

  % Check resizing total conserves power
  Ttotal = T;
  Ttotal.Nmax = Ttotal.Nmax + 5;
  sz = size(Ttotal.data, 2);
  testCase.assertEqual(Ttotal.columnCheck(), ones(1, sz), ...
    'AbsTol', 1.0e-15, 'Ttotal');

  % Check resizing scattered conserves power
  Tscat = T.scattered;
  Tscat.Nmax = Tscat.Nmax + 5;
  sz = size(Tscat.data, 2);
  testCase.assertEqual(Tscat.columnCheck(), ones(1, sz), ...
    'AbsTol', 1.0e-15, 'Tscat');
end

function testResizing(testCase)
  % Check resizing T-matrix works

  T = ott.tmatrix.Tmatrix(eye(6), 'type', 'scattered');
  Tnew1 = T;
  
  % Test same size
  Tnew1.Nmax = Tnew1.Nmax;
  testCase.verifyEqual(Tnew1, T, 'should not change');

  Tnew1.Nmax = Tnew1.Nmax + 5;
  testCase.verifyEqual(Tnew1.Nmax, T.Nmax + 5, ...
      'Faild to increase Nmax with vector size');
  testCase.verifyTrue(all(size(Tnew1.data) > size(T.data)), ...
      'Tmatrix size not actually increased (vector input)');

  Tnew2 = Tnew1;
  Tnew2.Nmax = T.Nmax;
  testCase.verifyEqual(Tnew2.Nmax, T.Nmax, ...
      'Failed to decrease Nmax with vector size');
  testCase.verifyEqual(size(Tnew2.data), size(T.data), ...
      'Tmatrix size not actually decreased (vector input)');

  Tnew1 = T;
  Tnew1.Nmax = T.Nmax(1) + 5;
  testCase.verifyEqual(Tnew1.Nmax, T.Nmax + 5, ...
      'Faild to increase Nmax (scalar input)');
  testCase.verifyTrue(all(size(Tnew1.data) > size(T.data)), ...
      'Tmatrix size not actually increased (scalar input)');

  Tnew1 = T;
  Tnew1.Nmax = [T.Nmax(1) + 5, T.Nmax(2)];
  testCase.verifyEqual(Tnew1.Nmax, [T.Nmax(1) + 5, T.Nmax(2)], ...
      'Faild to increase Nmax (uneven input)');
  testCase.verifyTrue(size(Tnew1.data, 1) > size(T.data, 1) ...
      && size(Tnew1.data, 2) == size(T.data, 2), ...
      'Tmatrix size not increased correctly (uneven input)');

  Tnew1 = T;
  Tnew1.Nmax(1) = T.Nmax(1) + 5;
  testCase.verifyEqual(Tnew1.Nmax, [T.Nmax(1) + 5, T.Nmax(2)], ...
      'Faild to increase Nmax (index input)');
  testCase.verifyTrue(size(Tnew1.data, 1) > size(T.data, 1) ...
      && size(Tnew1.data, 2) == size(T.data, 2), ...
      'Tmatrix size not increased correctly (index input)');
end

function testShrinkNmax(testCase)

  tmatrix = ott.tmatrix.Tmatrix(zeros(6));
  testCase.assertEqual(tmatrix.Nmax, [1, 1]);
  tmatrix = tmatrix.shrinkNmax();
  testCase.verifyEqual(tmatrix.Nmax, [0, 0]);

end

function testPlus(testCase)

  tmatrix1 = ott.tmatrix.Tmatrix(randn(6));
  tmatrix2 = ott.tmatrix.Tmatrix(randn(6));

  S = tmatrix1 + tmatrix2;
  testCase.verifyEqual(S.data, ...
    tmatrix1.data + tmatrix2.data, 'plus tmatrix');
  
  S = tmatrix1 + 1.0;
  testCase.verifyEqual(S.data, ...
    tmatrix1.data + 1.0, 'plus scalar');

end

function testMinus(testCase)

  tmatrix1 = ott.tmatrix.Tmatrix(randn(6));
  tmatrix2 = ott.tmatrix.Tmatrix(randn(6));

  S = tmatrix1 - tmatrix2;
  testCase.verifyEqual(S.data, ...
    tmatrix1.data - tmatrix2.data, 'minus tmatrix');
  
  S = 1 - tmatrix1;
  testCase.verifyEqual(S.data, ...
    1 - tmatrix1.data, 'minus scalar');

end

function testDivide(testCase)

  tmatrix = ott.tmatrix.Tmatrix(randn(6));

  S = tmatrix ./ 2;
  testCase.verifyEqual(S.data, tmatrix.data./2);

  S = tmatrix / 2;
  testCase.verifyEqual(S.data, tmatrix.data./2);

end

function testTime(testCase)

  tmatrix1 = ott.tmatrix.Tmatrix(randn(6));
  tmatrix2 = ott.tmatrix.Tmatrix(randn(6));
  S = tmatrix1 .* tmatrix2;
  testCase.verifyEqual(S.data, tmatrix1.data.*tmatrix2.data, 'times');
end

function testMtimes(testCase)

  tmatrix1 = ott.tmatrix.Tmatrix(randn(6));
  tmatrix2 = ott.tmatrix.Tmatrix(randn(6));
  
  S = tmatrix1 * tmatrix2;
  testCase.verifyEqual(S.data, tmatrix1.data * tmatrix2.data, 'T*T');
  
  beam = ott.bsc.Bsc(ones(3, 1), ones(3, 1));
  S = tmatrix1 * beam;
  testCase.verifyEqual(S.getCoefficients, sum(tmatrix1.data, 2), 'T*B');
  
  S = tmatrix1 * 2.0;
  testCase.verifyEqual(S.data, tmatrix1.data * 2.0, 'T*S');
  
  S = tmatrix1 .* 2.0;
  testCase.verifyEqual(S.data, tmatrix1.data * 2.0, 'T.*S');

end

function testDiag(testCase)

  dataA = randn(3);
  dataB = randn(3);
  dataC = randn(3);
  data = [dataA, dataC; -dataC, dataB];
  tmatrix = ott.tmatrix.Tmatrix(data);
  
  testCase.verifyEqual(diag(tmatrix), diag(data), 'diag');
  
  testCase.verifyEqual(diag(tmatrix, 1), ...
      [diag(dataA, 1); diag(dataB, 1)], 'diag k = 1');
  
  [da, db] = diag(tmatrix);
  testCase.verifySize(da, [3, 1], 'sz da');
  testCase.verifySize(db, [3, 1], 'sz db');
  testCase.verifyEqual([da; db], diag(data), '[da; db]');

end

function testMergeCols(testCase)

  tmatrix1 = ott.tmatrix.Tmatrix(zeros(6));
  tmatrix2 = ott.tmatrix.Tmatrix(eye(6));
  target = ott.tmatrix.Tmatrix(diag([0,0,1,0,0,1]));
  
  tmatrix = tmatrix1.mergeCols(tmatrix2, 3);
  testCase.verifyEqual(tmatrix, target);
end