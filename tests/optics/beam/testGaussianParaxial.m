function tests = testGaussianParaxial
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  waist = 1.0;
  beam = ott.optics.beam.GaussianParaxial(waist);
  
end

function testPoynting(testCase)

  waist = 1.0;
  beam = ott.optics.beam.GaussianParaxial(waist);
  
  xyz = randn(3, 5);
  S = beam.poynting(xyz);
  
  testCase.verifyEqual(S.vxyz(1:2, :), zeros(2, 5), 'zero paraxial');
  testCase.verifyNotEqual(S.vxyz(3, :), zeros(1, 5), 'z');

end

function testVisualise(testCase)

  waist = 1.0;
  beam = ott.optics.beam.GaussianParaxial(waist);
  
  h = figure();
  beam.visualise();
  close(h);
  
end
