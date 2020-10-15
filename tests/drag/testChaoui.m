function tests = testFaxen
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');

  % Turn off warning about FaxnTerms
  warning('off', 'ott:drag:ChaouiSphere:faxen_perp_terms');

  % Turn on intended warnings
  warning('on', 'ott:drag:ChaouiSphere:large_epsilon');
end

function testConstruction(testCase)

  radius = 1.0;
  viscosity = 1.4;
  separation = 1.2;
  a = ott.drag.ChaouiSphere(radius, separation, viscosity);

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

function testLargeEpsilonWarning(testCase)

  radius = 1.0;
  viscosity = 1.0;
  separation = 1e6;
  drag = ott.drag.ChaouiSphere(radius, separation, viscosity);

  testCase.verifyWarning(@() drag.forward, ...
    'ott:drag:ChaouiSphere:large_epsilon');

end

function testFarLimit(testCase)

  radius = 1.0;
  viscosity = 1.0;
  separation = radius + 0.2*radius;

  a = ott.drag.ChaouiSphere(radius, separation, viscosity);
  b = ott.drag.PadeSphere(radius, separation, viscosity);

  S = warning('off', 'ott:drag:ChaouiSphere:large_epsilon');

  testCase.verifyEqual(a.forward, b.forward, ...
      'RelTol', 6.0e-2, ...
      'Far-limit doesn''t match forward');

  warning(S);

end

function testPlot(testCase)
  % Attempt to reproduce figure 1 from M. Chaoui and F. Feuillebois
  % https://doi.org/10.1093/qjmam/56.3.381

  radius = 1.0;
  separation = radius + radius.*logspace(-6, 1, 100);
  viscosity = 1.0;

  values = zeros(3, length(separation));

  sphere = ott.drag.StokesSphere(radius, viscosity);

  S1 = warning('off', 'ott:drag:ChaouiSphere:large_epsilon');
  S2 = warning('off', 'ott:drag:PadeSphere:small_epsilon');
  S3 = warning('off', 'ott:drag:FaxenSphere:small_epsilon');

  idx = sub2ind(size(sphere.forward), 1, 1);
%   idx = sub2ind(size(sphere.forward), 5, 1);  % NO!
%   idx = sub2ind(size(sphere.forward), 4, 4);

  for ii = 1:length(separation)
    drag = ott.drag.ChaouiSphere(radius, separation(ii), viscosity);
    values(1, ii) = drag.forward(idx);

    drag = ott.drag.FaxenSphere(radius, separation(ii), viscosity);
    values(2, ii) = drag.forward(idx);

    drag = ott.drag.PadeSphere(radius, separation(ii), viscosity);
    values(3, ii) = drag.forward(idx);
  end

  warning(S1);
  warning(S2);
  warning(S3);

%   figure();
%   semilogx(separation./radius - 1, values./sphere.forward(idx));
%   axis([1e-6, 1e1, 1, 8]);

%   figure();
%   loglog(separation./radius - 1, abs((values(2, :) - values(1, :))./values(1, :)));

end

