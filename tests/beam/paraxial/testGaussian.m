function tests = testGaussianParaxial
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  waist = 1.0;
  beam = ott.beam.paraxial.Gaussian(waist);
  
  testCase.verifyEqual(beam.waist, waist, 'waist');
  testCase.verifyEqual(beam.power, 1.0, 'power');
  testCase.verifyEqual(beam.polarisation, [1, 0], 'pol');
  
end

function testAgainstDavis5(testCase)

  waist = 5.0;
  beam0 = ott.beam.paraxial.Gaussian(waist);
  beam5 = ott.beam.GaussianDavis5(waist);
  
  xyz = randn(3, 5);
  [E0, H0] = beam0.ehfield(xyz);
  [E5, H5] = beam5.ehfield(xyz);
  
  testCase.verifyEqual(E0.vxyz, E5.vxyz, ...
    'AbsTol', 2e-2, 'E fields don''t match');
  testCase.verifyEqual(H0.vxyz, H5.vxyz, ...
    'AbsTol', 2e-2, 'H fields don''t match');
end

function testPoynting(testCase)

  waist = 1.0;
  beam = ott.beam.paraxial.Gaussian(waist);
  
  xyz = randn(3, 5);
  S = beam.poynting(xyz);
  
  testCase.verifyEqual(S.vxyz(1:2, :), zeros(2, 5), 'zero paraxial');
  testCase.verifyNotEqual(S.vxyz(3, :), zeros(1, 5), 'z');

end

function testVisualise(testCase)

  waist = 1.0;
  beam = ott.beam.paraxial.Gaussian(waist);
  
  h = figure();
  beam.visualise();
  close(h);
  
end

function testVisualiseFarfieldSphere(testCase)

  waist = 0.5;
  beam = ott.beam.paraxial.Gaussian(waist);
  
  h = figure();
  beam.visualiseFarfieldSphere();
  close(h);
  
end
