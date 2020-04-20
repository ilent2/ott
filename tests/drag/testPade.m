function tests = testPade
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruction(testCase)

  radius = 1.0;
  viscosity = 1.4;
  separation = 1.2;
  a = ott.drag.PadeSphere(radius, separation, viscosity);

  testCase.verifyEqual(a.separation, separation, 'Separation not set');
  testCase.verifyEqual(a.radius, radius, 'Separation not set');
  testCase.verifyEqual(a.viscosity, viscosity, 'Separation not set');

  testCase.verifyNotEqual(a.forward(1, 1), a.forward(3, 3), ...
    'A(1,1) should not equal A(3,3)');
  testCase.verifyEqual(a.forward(1, 1), a.forward(2, 2), ...
    'A(1,1) should equal A(2,2)');
  testCase.verifyEqual(a.forward(4, 4), a.forward(5, 5), ...
    'A(4,4) should equal A(5,5)');
  testCase.verifyNotEqual(a.forward(5, 5), a.forward(6, 6), ...
    'A(5,5) should equal A(6,6)');
  
  testCase.verifyNotEqual(a.forward(1, 5), 0, ...
    'A(1, 5) should not be zero');
  testCase.verifyNotEqual(a.forward(2, 4), 0, ...
    'A(2, 4) should not be zero');
  testCase.verifyNotEqual(a.forward(4, 2), 0, ...
    'A(4, 2) should not be zero');
  testCase.verifyNotEqual(a.forward(5, 1), 0, ...
    'A(5, 1) should not be zero');

end

function testSmallEpsilonWarning(testCase)

  radius = 1.0;
  viscosity = 1.0;
  separation = radius .* (1 + 1e-4);
  
  testCase.verifyWarning(@() ott.drag.PadeSphere(radius, separation, viscosity), ...
    'ott:drag:PadeSphere:small_epsilon');

end

function testFarLimit(testCase)

  radius = 1.0;
  viscosity = 1.0;
  separation = radius .* (1 + 1);
  
  S = warning('off', 'ott:drag:PadeSphere:small_epsilon');

  a = ott.drag.PadeSphere(radius, separation, viscosity);
  b = ott.drag.FaxenSphere(radius, separation, viscosity);
  
  warning(S);

  testCase.verifyEqual(a.forward, b.forward, ...
      'RelTol', 6.0e-2, ...
      'Far-limit doesn''t match forward');
  testCase.verifyEqual(a.inverse, b.inverse, ...
      'RelTol', 5.0e-2, ...
      'Far-limit doesn''t match inverse');

end

