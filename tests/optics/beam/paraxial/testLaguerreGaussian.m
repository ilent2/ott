function tests = testLaguerreGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../../');
end

function testConstructor(testCase)

  waist = 1.0;
  lmode = 3;
  pmode = 2;
  beam = ott.optics.beam.paraxial.LaguerreGaussian(waist, lmode, pmode);
  testCase.verifyEqual(beam.waist, waist, 'w');
  testCase.verifyEqual(beam.lmode, lmode, 'l');
  testCase.verifyEqual(beam.pmode, pmode, 'p');
end

function testZerothOrder(testCase)

  waist = 1.0;
  beam0 = ott.optics.beam.paraxial.LaguerreGaussian(waist, 0, 0);
  beamG = ott.optics.beam.paraxial.Gaussian(waist);
  
  xyz = [randn(2, 5); zeros(1, 5)];
%   xyz = randn(3, 5);
  E0 = beam0.efield(xyz);
  EG = beamG.efield(xyz);
  
  testCase.verifyEqual(E0.vxyz, EG.vxyz, ...
    'AbsTol', 1e-15, 'fields don''t match');
end

function testVisualise(testCase)

  waist = 1.0;
  lmode = 2;
  pmode = 3;
  beam = ott.optics.beam.paraxial.LaguerreGaussian(waist, lmode, pmode);
  
  h = figure();
  beam.visualise('range', [3, 3].*beam.waist);
  close(h);

end
