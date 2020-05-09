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
 
  index_relative = 1.0;
  plane = ott.scat.planewave.Plane([-1;0;0], index_relative);
  beam = ott.beam.abstract.PlaneWave('direction', [1;0;0], ...
    'polarisation', [0;1;0], 'origin', [0;0;0], ...
    'field', 1.0);
  
  [rbeam, tbeam] = plane.scatter(beam);
  
  testCase.verifyEqual(tbeam.field, [0;1].*ones(size(beam, 2)), ...
    'tbeam = ibeam field');
  testCase.verifyEqual(tbeam.direction, beam.direction, ...
    'tbeam = ibeam direction');
  testCase.verifyEqual(rbeam.field, [0.0; 0.0]);

end

function testScatterDielectric(testCase)
  
  % Compare against fresnel coefficients (check vector math in scatter)
  index_relative = 1.5;
  wavenumber = 2*pi;
  kix = wavenumber;
  ktx = index_relative*wavenumber;
  [Sr, St] = ott.scat.planewave.Plane.fresnelS(kix, ktx, 1.0, index_relative);
 
  beam = ott.beam.abstract.PlaneWave('direction', [1;0;0], ...
    'polarisation', [0;1;0], 'origin', [0;0;0], ...
    'field', 1.0, 'wavenumber', wavenumber);
  plane = ott.scat.planewave.Plane([-1;0;0], index_relative);
  [rbeam, tbeam] = plane.scatter(beam);
  
  testCase.verifyEqual(rbeam.field, [0;Sr].*beam.field, 'rbeam');
  testCase.verifyEqual(tbeam.field, [0;St].*beam.field, 'tbeam');
  testCase.verifyEqual(rbeam.wavenumber, 2*pi, 'orel');
  testCase.verifyEqual(tbeam.wavenumber, 2*pi*index_relative, 'irel');
 
  beam = ott.beam.abstract.PlaneWave('direction', [-1;0;0], ...
    'polarisation', [0;1;0], 'origin', [0;0;0], ...
    'field', 1.0, 'wavenumber', wavenumber);
  plane = ott.scat.planewave.Plane([1;0;0], index_relative);
  [rbeam, tbeam] = plane.scatter(beam);
  
  testCase.verifyEqual(rbeam.field, [0;Sr].*beam.field, 'rbeam');
  testCase.verifyEqual(tbeam.field, [0;St].*beam.field, 'tbeam');

end

function testScatterReflection(testCase)
  
  % Compare against fresnel coefficients (check vector math in scatter)
  index_relative = 1.5;
  wavenumber = 2*pi;
  kix = index_relative*wavenumber;
  ktx = wavenumber;
  [Sr, St] = ott.scat.planewave.Plane.fresnelS(kix, ktx, index_relative, 1.0);
 
  beam = ott.beam.abstract.PlaneWave('direction', [1;0;0], ...
    'polarisation', [0;1;0], 'origin', [0;0;0], ...
    'field', 1.0, 'index', index_relative, 'wavenumber', index_relative*wavenumber);
  plane = ott.scat.planewave.Plane([1;0;0], index_relative);
  [rbeam, tbeam] = plane.scatter(beam);
  
  testCase.verifyEqual(rbeam.field, [0;Sr].*beam.field, ...
    'AbsTol', 1.0e-15, 'rbeam');
  testCase.verifyEqual(tbeam.field, [0;St].*beam.field, ...
    'AbsTol', 1.0e-15, 'tbeam');
  testCase.verifyEqual(rbeam.medium.index, index_relative, 'irel');
  testCase.verifyEqual(tbeam.medium.index, 1.0, 'orel');

end

function testPowerChange(testCase)

  index_relative = 1.5;
  wavenumber = 2*pi;
  kix = wavenumber;
  ktx = index_relative*wavenumber;
  [Sr, St] = ott.scat.planewave.Plane.fresnelS(kix, ktx, 1.0, index_relative);
 
  beam = ott.beam.abstract.PlaneWave('direction', [1;0;0], ...
    'polarisation', [0;1;0], 'origin', [0;0;0], ...
    'field', 0.5, 'wavenumber', wavenumber);
  plane = ott.scat.planewave.Plane([-1;0;0], index_relative);
  [rbeam, tbeam] = plane.scatter(beam);
  
  testCase.verifyEqual(rbeam.field, [0;Sr].*beam.field, 'rbeam');
  testCase.verifyEqual(tbeam.field, [0;St].*beam.field, 'tbeam');
end

function testArrayOfArrayWaves(testCase)

  index_relative = 1.5;

  beam1 = ott.beam.abstract.PlaneWave('direction', [1;0;0], ...
    'polarisation', [0;1;0], 'origin', [0;0;0], ...
    'field', 0.5, 'index', 1.0);
  beam2 = ott.beam.abstract.PlaneWave('direction', [-1;0;0], ...
    'polarisation', [0;1;0], 'origin', [0;0;0], ...
    'field', 0.5, 'index', index_relative);
  beam = [beam1, beam2];
  
  plane = ott.scat.planewave.Plane([-1;0;0], index_relative);
  
  [rbeam, tbeam] = plane.scatter(beam);
  
  testCase.verifyClass(rbeam, 'ott.beam.Array');
  testCase.verifyClass(tbeam, 'ott.beam.Array');
  testCase.verifyEqual(size(rbeam), [1, 2], 'rbeam sz');
  testCase.verifyEqual(size(tbeam), [1, 2], 'tbeam sz');
  testCase.verifyEqual(rbeam(1).medium.index, 1.0, 'rbeam1 idx');
  testCase.verifyEqual(rbeam(2).medium.index, index_relative, 'rbeam2 idx');
  testCase.verifyEqual(tbeam(1).medium.index, index_relative, 'tbeam1 idx');
  testCase.verifyEqual(tbeam(2).medium.index, 1.0, 'tbeam2 idx');
end