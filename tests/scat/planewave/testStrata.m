function tests = testStrata
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testCoefficientsIndexMatched(testCase)

  n1 = 1;
  n2 = 1;
  kix = n1*2*pi;
  ktx = n2*2*pi;
  
  depths = [];
  normal = [0;0;1];
  direction = [0;1;1];
  kvecs = [kix, ktx] .* direction ./ vecnorm(direction);
  
  % Calculate S and P coefficients
  [aS, bS] = ott.scat.planewave.Strata.coeffS(normal, depths, kvecs);
  [aP, bP] = ott.scat.planewave.Strata.coeffP(normal, depths, kvecs);
  
  testCase.verifyEqual(bS, 0.0, 'Sr n1 = n2');
  testCase.verifyEqual(aS, 1.0, 'St n1 = n2');
  
  testCase.verifyEqual(bP, 0.0, 'Pr n1 = n2');
  testCase.verifyEqual(aP, 1.0, 'Pt n1 = n2');

end

function testCoefficientsIndexMatched3Layeres(testCase)

  n1 = 1;
  n2 = 1;
  n3 = 1;
  ks = [n1, n2, n3]*2*pi;
  
  depths = [0.2];     % [L]
  normal = [0;0;1];
  direction = [0;1;1];
  kvecs = ks .* direction ./ vecnorm(direction);    % [2*pi/L]
  
  % Calculate S and P coefficients
  [aS, bS] = ott.scat.planewave.Strata.coeffS(normal, depths, kvecs);
  [aP, bP] = ott.scat.planewave.Strata.coeffP(normal, depths, kvecs);
  
  testCase.verifyEqual(bS, [0.0; 0.0], 'AbsTol', 1.0e-14, 'Sr n1 = n2');
  testCase.verifyEqual(aS, [1.0; 1.0], 'AbsTol', 1.0e-14, 'St n1 = n2');
  
  testCase.verifyEqual(bP, [0.0; 0.0], 'AbsTol', 1.0e-14, 'Pr n1 = n2');
  testCase.verifyEqual(aP, [1.0; 1.0], 'AbsTol', 1.0e-14, 'Pt n1 = n2');

end

function testCompareCoefficientsToFresnel(testCase)
  % Compare to Plane.fresnel with a single interface

  n1 = 1;
  n2 = 1.3;
  kix = n1*2*pi;
  
  depths = [];
  normal = -[0;0;1];
  direction = [0;0.5;1];
  ki = kix .* direction ./ vecnorm(direction);
  kt = [ki(1:2); sqrt(vecnorm(ki).^2 .* (n2 ./ n1).^2 - vecnorm(ki(1:2)).^2)];
  kvecs = [ki, kt];
  
  % Calcualte Fresnel coefficients
  [Sr, St] = ott.scat.planewave.Plane.fresnelS(kvecs(3, 1), kvecs(3, 2), n1, n2);
  [Pr, Pt] = ott.scat.planewave.Plane.fresnelP(kvecs(3, 1), kvecs(3, 2), n1, n2);
  
  % Calculate S and P coefficients
  [aS, bS] = ott.scat.planewave.Strata.coeffS(normal, depths, kvecs);
  [aP, bP] = ott.scat.planewave.Strata.coeffP(normal, depths, kvecs);
  
  testCase.verifyEqual(bS, Sr, 'AbsTol', 1.0e-15, 'Sr Fres');
  testCase.verifyEqual(aS, St, 'AbsTol', 1.0e-15, 'St Fres');
  
  testCase.verifyEqual(bP, Pr, 'AbsTol', 1.0e-15, 'Pr Fres');
  testCase.verifyEqual(aP, Pt, 'AbsTol', 1.0e-15, 'Pt Fres');

end

function testCoefficientPureComplex(testCase)

  % Normal incidence, medium 2 is perfectly absorbing
  n1 = 1;
  n2 = -1i;
  kix = n1*2*pi;
  ktx = n2*2*pi;
  
  depths = [];
  normal = [0;0;1];
  direction = normal;
  kvecs = [kix, ktx] .* direction ./ vecnorm(direction);
  
  % Calculate S and P coefficients
  [aS, bS] = ott.scat.planewave.Strata.coeffS(normal, depths, kvecs);
  [aP, bP] = ott.scat.planewave.Strata.coeffP(normal, depths, kvecs);
  
  testCase.verifyEqual(abs(bS).^2, 1.0, 'AbsTol', 1.0e-15, 'Sr n2 = -1i');
  testCase.verifyTrue(~isreal(aS), 'St n2 = -1i');
  
  testCase.verifyEqual(abs(bP).^2, 1.0, 'AbsTol', 1.0e-15, 'Pr n2 = -1i');
  testCase.verifyTrue(~isreal(aP), 'Pt n2 = -1i');

end

