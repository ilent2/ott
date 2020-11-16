function tests = testFieldVector
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testCartesianToSpherical(testCase)

  xyz = randn(3, 5);
  xyzv = randn(3, 5);
  S = ott.utils.FieldVectorCart(xyzv, xyz);

  [rtpv, rtp] = ott.utils.xyzv2rtpv(xyzv, xyz);
  
  % Check original
  testCase.verifyEqual(S.vxyz, xyzv, 'xyzv');
  testCase.verifyEqual(S.vrtp, rtpv, 'rtpv');
  testCase.verifyEqual(double(S(4:6, :)), xyz, 'xyz');
  
  % Cast to Sph
  P = ott.utils.FieldVectorSph(S);
  testCase.verifyEqual(P.vxyz, xyzv, 'AbsTol', 1e-13, 'xyzv cast');
  testCase.verifyEqual(P.vrtp, rtpv, 'AbsTol', 1e-13, 'rtpv cast');
  testCase.verifyEqual(double(P(4:6, :)), rtp, 'rtp');
  
end

function testSphericalToCartesian(testCase)

  rtp = [abs(randn(1, 5)); pi*rand(1, 5); 2*pi*rand(1, 5)];
  rtpv = randn(3, 5);
  S = ott.utils.FieldVectorSph(rtpv, rtp);

  [xyzv, xyz] = ott.utils.rtpv2xyzv(rtpv, rtp);
  
  % Check original
  testCase.verifyEqual(S.vrtp, rtpv, 'rtpv');
  testCase.verifyEqual(S.vxyz, xyzv, 'xyzv');
  testCase.verifyEqual(double(S(4:6, :)), rtp, 'rtp');
  
  % Cast to Cart
  P = ott.utils.FieldVectorCart(S);
  testCase.verifyEqual(P.vxyz, xyzv, 'AbsTol', 1e-13, 'xyzv cast');
  testCase.verifyEqual(P.vrtp, rtpv, 'AbsTol', 1e-13, 'rtpv cast');
  testCase.verifyEqual(double(P(4:6, :)), xyz, 'xyz');
  
end

function testSphOutsideRange(testCase)

  rtp = [[-1;0.2;0.3], [1;-0.5;0.0], [1;0.2;5*pi]];
  rtpv = randn(3, 3);
  
  testCase.verifyWarning(@construct, ...
    'ott:utils:FieldVectorSph:rtp_outside_range');
  
  function S = construct
    S = ott.utils.FieldVectorSph(rtpv, rtp);
  end

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

