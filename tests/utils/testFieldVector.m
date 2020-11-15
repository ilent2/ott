function tests = testFieldVector
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testCartesianToSpherical(testCase)

  loc = randn(3, 5);
  val = randn(3, 5);
  S = ott.utils.FieldVectorCart(val, loc);

  rtp = ott.utils.xyzv2rtpv(val, loc);
  
  testCase.verifyEqual(S.vxyz, val, 'xyz');
  testCase.verifyEqual(S.vrtp, rtp, 'xyz');
  
end

function testSphericalToCartesian(testCase)

  loc = randn(3, 5);
  val = randn(3, 5);
  S = ott.utils.FieldVectorSph(val, loc);

  xyz = ott.utils.rtpv2xyzv(val, loc);
  
  testCase.verifyEqual(S.vrtp, val, 'xyz');
  testCase.verifyEqual(S.vxyz, xyz, 'xyz');
  
end

function testMathsCart(testCase)

  v1 = [1;2;3];
  v2 = [3;4;5];
  fv1 = ott.utils.FieldVectorCart(v1);
  fv2 = ott.utils.FieldVectorCart(v2);
  
  s1 = fv1 + fv2;
  testCase.verifyInstanceOf(s1, 'ott.utils.FieldVectorCart', '++');
  testCase.verifyEqual(double(s1), v1 + v2, '++');
  
  s1 = 3 + fv2;
  testCase.verifyInstanceOf(s1, 'double', '++3');
  testCase.verifyEqual(double(s1), 3 + v2, '++3');
  
  s1 = fv1 + 5;
  testCase.verifyInstanceOf(s1, 'double', '++5');
  testCase.verifyEqual(double(s1), v1 + 5, '++5');
  
  s1 = fv1 - fv2;
  testCase.verifyInstanceOf(s1, 'ott.utils.FieldVectorCart', '--');
  testCase.verifyEqual(double(s1), v1 - v2, '--');
  
  s1 = 3 - fv2;
  testCase.verifyInstanceOf(s1, 'double', '--3');
  testCase.verifyEqual(double(s1), 3 - v2, '--3');
  
  s1 = fv1 - 5;
  testCase.verifyInstanceOf(s1, 'double', '--5');
  testCase.verifyEqual(double(s1), v1 - 5, '--5');
  
  s1 = - fv2;
  testCase.verifyInstanceOf(s1, 'ott.utils.FieldVectorCart', '-');
  testCase.verifyEqual(double(s1), - v2, 'sum');
  
  s1 = 1.5 * fv2;
  testCase.verifyInstanceOf(s1, 'ott.utils.FieldVectorCart', '*');
  testCase.verifyEqual(double(s1), 1.5 * v2, 'sum');
  
  s1 = 1.5 .* fv2;
  testCase.verifyInstanceOf(s1, 'ott.utils.FieldVectorCart', '.*');
  testCase.verifyEqual(double(s1), 1.5 * v2, 'sum');
  
  s1 = fv1 / 2;
  testCase.verifyInstanceOf(s1, 'ott.utils.FieldVectorCart', '/');
  testCase.verifyEqual(double(s1), v1 / 2, 'sum');
  
  s1 = fv1 ./ 2;
  testCase.verifyInstanceOf(s1, 'ott.utils.FieldVectorCart', './');
  testCase.verifyEqual(double(s1), v1 / 2, 'sum');

end

function testSubsasgnDifferentTypes(testCase)

  loc = randn(3, 6);
  val = randn(3, 6);
  Cart = ott.utils.FieldVectorCart(val, loc);
  Sph = ott.utils.FieldVectorSph(Cart);
  
  % Explicit cast
  cpy = Cart;
  cpy(:, 4:6) = ott.utils.FieldVectorCart(Sph(:, 4:6));
  testCase.verifyEqual(double(cpy), double(Cart), 'AbsTol', 1e-14, ...
    'explicit cast');
  
  % Explicit cast with 'like'
  % Doesn't seem to work by default in 2018a or 2020b, added a workaround
  % to FieldVector class (cast method) to implement the behaviour I expect
  cpy = Cart;
  cpy(:, 4:6) = cast(Sph(:, 4:6), 'like', cpy);
  testCase.verifyEqual(double(cpy), double(Cart), 'AbsTol', 1e-14, ...
    'explicit cast with like');
  
  % Implicit cast
  % Doesn't seem to work in 2018a, works fine in 2020b
  cpy = Cart;
  cpy(:, 4:6) = Sph(:, 4:6);
  testCase.verifyEqual(double(cpy), double(Cart), 'AbsTol', 1e-14, ...
    'implicit cast');

end

