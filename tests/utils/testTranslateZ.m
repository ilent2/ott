function tests = translate_z
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testValues(testCase)
  % Test values in the matricies

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-3;

  % This should return sparse eye matrics
  [Aeye, Beye] = ott.utils.translate_z(5, 0.0);

  sz = ott.utils.combined_index(5, 5);
  testCase.verifyThat(size(Aeye), IsEqualTo([sz, sz]), ...
      'Incorrect matrix size');
  testCase.verifyThat(Aeye, IsEqualTo(speye([sz, sz]), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect A matrix values');
  testCase.verifyThat(Beye, IsEqualTo(sparse(sz, sz), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect B matrix values');

  % Check some values
  [A, B] = ott.utils.translate_z(5, pi/2);

  targetA = sparse(1, [1, 5, 11, 19, 29], ...
      [-0.0786, 0.1411, 0.1964, -0.0217, -0.2398], 1, sz);
  testCase.verifyThat(A(1, :), IsEqualTo(targetA, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect A matrix values');
  targetB = sparse(1, [1, 5, 11, 19, 29], ...
      [-0.1306i, -0.1357i, 0.1181i, 0.2770i, 0.1312i], 1, sz);
  testCase.verifyThat(B(1, :), IsEqualTo(targetB, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect B matrix values');

end

function testDifferentNmax(testCase)
  % Test outputs with different row/column Nmax

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-3;

  [Afull, Bfull] = ott.utils.translate_z(7, pi/2);

  sz1 = ott.utils.combined_index(5, 5);
  sz2 = ott.utils.combined_index(7, 7);

  [A, B] = ott.utils.translate_z([5, 7], pi/2);
  testCase.verifyThat(size(A), IsEqualTo([sz1, sz2]), ...
      'Incorrect matrix size [5 7]');
  testCase.verifyThat(A, IsEqualTo(Afull(1:sz1, 1:sz2), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect A matrix values [5 7]');
  testCase.verifyThat(B, IsEqualTo(Bfull(1:sz1, 1:sz2), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect B matrix values [5 7]');

  [A, B] = ott.utils.translate_z([7, 5], pi/2);
  testCase.verifyThat(size(A), IsEqualTo([sz2, sz1]), ...
      'Incorrect matrix size [7 5]');
  testCase.verifyThat(A, IsEqualTo(Afull(1:sz2, 1:sz1), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect A matrix values [7 5]');
  testCase.verifyThat(B, IsEqualTo(Bfull(1:sz2, 1:sz1), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect B matrix values [7 5]');
end

function testLargeTranslations(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  % For this translation the beam power should go to zero

  beam = ott.BscPmGauss();
  beam.power = 1.0;
  tbeam = beam.translateXyz([300;0;0]);  % calls translate_z

  testCase.verifyThat(tbeam.power, IsEqualTo(0.0, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Beam power does not drop to zero for large radial translations');

end

function testNegativeTranslations(testCase)
  
  Nmax = 1;
  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1e-5;
  
  % Test Gumerov
  
  [A1, B1] = ott.utils.translate_z(Nmax, 1e-6, 'method', 'gumerov');
  [A2, B2] = ott.utils.translate_z(Nmax, -1e-6, 'method', 'gumerov');
  C1 = [ A1 B1 ; B1 A1 ];
  C2 = [ A2 B2 ; B2 A2 ];
  
  testCase.verifyThat(C1*C2, IsEqualTo(speye(size(C1*C2)), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Gumerov fails negative translation');
  
  % Test Videen
  
  [A1, B1] = ott.utils.translate_z(Nmax, 1e-6, 'method', 'videen');
  [A2, B2] = ott.utils.translate_z(Nmax, -1e-6, 'method', 'videen');
  C1 = [ A1 B1 ; B1 A1 ];
  C2 = [ A2 B2 ; B2 A2 ];
  
  testCase.verifyThat(C1*C2, IsEqualTo(speye(size(C1*C2)), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Videen fails negative translation');
    
end

function testNegativeDifferentNmax(testCase)
  % Test outputs with different row/column Nmax

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-3;
  
  dz = -pi/2;

  [Afull, Bfull] = ott.utils.translate_z(7, dz);

  sz1 = ott.utils.combined_index(5, 5);
  sz2 = ott.utils.combined_index(7, 7);

  [A, B] = ott.utils.translate_z([5, 7], dz);
  testCase.verifyThat(size(A), IsEqualTo([sz1, sz2]), ...
      'Incorrect matrix size [5 7]');
  testCase.verifyThat(A, IsEqualTo(Afull(1:sz1, 1:sz2), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect A matrix values [5 7]');
  testCase.verifyThat(B, IsEqualTo(Bfull(1:sz1, 1:sz2), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect B matrix values [5 7]');

  [A, B] = ott.utils.translate_z([7, 5], dz);
  testCase.verifyThat(size(A), IsEqualTo([sz2, sz1]), ...
      'Incorrect matrix size [7 5]');
  testCase.verifyThat(A, IsEqualTo(Afull(1:sz2, 1:sz1), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect A matrix values [7 5]');
  testCase.verifyThat(B, IsEqualTo(Bfull(1:sz2, 1:sz1), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect B matrix values [7 5]');
end

