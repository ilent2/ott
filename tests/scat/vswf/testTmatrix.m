function tests = testTmatrix
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructDefault(testCase)

  T = ott.scat.vswf.Tmatrix();
  
  testCase.verifyEmpty(T.data, 'empty');
  testCase.verifyEqual(T.type, 'scattered', 'type');
  testCase.verifyEqual(T.position, [0;0;0], 'position');
  testCase.verifyEqual(T.rotation, eye(3), 'rotation');
end

function testEmpty(testCase)

  T = ott.scat.vswf.Tmatrix.empty();
  testCase.verifyEmpty(T);
end

function testScatteredToTotal(testCase)

  tmatrix = ott.scat.vswf.Tmatrix([], 'type', 'scattered');
  testCase.assertEqual(tmatrix.Nmax, [0,0], 'original Nmax');
  testCase.assertEqual(tmatrix.type, 'scattered');
  
  total = tmatrix.total;
  testCase.verifyEqual(total.data, [], 'empty');
  testCase.verifyEqual(total.type, 'total', 'new type');
  
  total.Nmax = 1;
  testCase.verifyEqual(total.data, eye(6), 'Nmax = 1');

end

function testTotalToScattered(testCase)

  tmatrix = ott.scat.vswf.Tmatrix(eye(6), 'type', 'total');
  testCase.assertEqual(tmatrix.Nmax, [1,1], 'original Nmax');
  testCase.assertEqual(tmatrix.type, 'total');
  
  scattered = tmatrix.scattered;
  testCase.verifyEqual(scattered.Nmax, [1,1], 'Nmax');
  testCase.verifyEqual(scattered.data, zeros(6), 'data');
  testCase.verifyEqual(scattered.type, 'scattered', 'new type');
  
  scattered.Nmax = 0;
  testCase.verifyEqual(scattered.data, [], 'empty');

end

function testRealImagFunctions(testCase)

  data = randn(6, 6) + 1i*randn(6, 6);
  T = ott.scat.vswf.Tmatrix(data);

  Treal = real(T);
  testCase.verifyEqual(Treal.data, real(data), 'real');
  
  Timag = imag(T);
  testCase.verifyEqual(Timag.data, imag(data), 'real');

end

function testShrinkPowerWarning(testCase)
  % Check shrinking gives a power warning

  T = ott.scat.vswf.Tmatrix(eye(16), 'type', 'scattered');
  
  testCase.verifyWarning(@() T.setNmax(1), ...
      'ott:Tmatrix:setNmax:truncation', 'warn');
    
  testCase.verifyError(@() T.setNmax(1, 'powerloss', 'error'), ...
      'ott:Tmatrix:setNmax:truncation', 'error');

  testCase.verifyWarningFree(@() T.setNmax(1, 'powerloss', 'ignore'), ...
      'ignore');
end

function testPowerOnResize(testCase)
  
  % Check columnCheck
  T = ott.scat.vswf.Tmatrix(eye(6), 'type', 'total');
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

  T = ott.scat.vswf.Tmatrix(eye(6), 'type', 'scattered');
  Tnew1 = T;

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
