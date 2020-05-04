function tests = testPlane
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testFresnelCoefficientsIndexMatched(testCase)

  n1 = 1;
  n2 = 1;
  kix = n1*2*pi;
  ktx = n2*2*pi;
  
  [Sr, St] = ott.scat.planewave.Plane.fresnelS(kix, ktx, n1, n2);
  testCase.verifyEqual(Sr, 0.0, 'Sr n1 = n2');
  testCase.verifyEqual(St, 1.0, 'St n1 = n2');
  
  [Pr, Pt] = ott.scat.planewave.Plane.fresnelP(kix, ktx, n1, n2);
  testCase.verifyEqual(Pr, 0.0, 'Pr n1 = n2');
  testCase.verifyEqual(Pt, 1.0, 'Pt n1 = n2');

end

function testFresnelCoefficientsPureComplex(testCase)

  % Normal incidence, medium 2 is perfectly absorbing
  n1 = 1;
  n2 = -1i;
  kix = n1*2*pi;
  ktx = n2*2*pi;
  
  [Sr, St] = ott.scat.planewave.Plane.fresnelS(kix, ktx, n1, n2);
  testCase.verifyEqual(abs(Sr).^2, 1.0, 'Sr n2 = -1i');
  testCase.verifyTrue(~isreal(St), 'St n2 = -1i');
  
  [Pr, Pt] = ott.scat.planewave.Plane.fresnelP(kix, ktx, n1, n2);
  testCase.verifyEqual(abs(Pr).^2, 1.0, 'Pr n2 = -1i');
  testCase.verifyTrue(~isreal(Pt), 'Pt n2 = -1i');

end

function testFresnelCoefficientsTir(testCase)

  % 45 deg incidence, medium 2 is much higher refractive index
  n1 = 1;
  n2 = 2;
  ang = pi/4;
  kix = 2*pi*sin(ang);
  ktx = 2*pi*sqrt(1 - (n2./n1).^2 .* sin(ang).^2);
  
  [Sr, St] = ott.scat.planewave.Plane.fresnelS(kix, ktx, n1, n2);
  testCase.verifyEqual(abs(Sr).^2, 1.0, 'AbsTol', 1e-15, 'Sr n2 = 2');
  testCase.verifyTrue(~isreal(St), 'St n2 = 2');
  
  [Pr, Pt] = ott.scat.planewave.Plane.fresnelP(kix, ktx, n1, n2);
  testCase.verifyEqual(abs(Pr).^2, 1.0, 'Pr n2 = 2');
  testCase.verifyTrue(~isreal(Pt), 'Pt n2 = 2');

end

function testScatterIndexMatched(testCase)

  plane = ott.scat.planewave.Plane(1.0, [1;0;0]);
  beam = ott.beam.abstract.PlaneWave('direction', [1;0;0], ...
    'polarisation', [0;1;0], 'origin', [0;0;0], ...
    'field', 1.0);
  
  [rbeam, tbeam] = plane.scatter(beam);
  
  testCase.verifyEqual(tbeam, beam, 'tbeam = ibeam');
  testCase.verifyEqual(rbeam.field, 0.0);

end
