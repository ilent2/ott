function tests = testGaussianDavis5
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstructor(testCase)

  waist = 1.0;
  beam = ott.beam.GaussianDavis5(waist);
  
end

function testPosition(testCase)

  waist = 1.0;
  beam = ott.beam.GaussianDavis5(waist);
  
  xyz = [0.1;0.3;-1];
  E0 = beam.efield([0;0;0]);
  
  beam.position = xyz;
  E1 = beam.efield(xyz);
  
  testCase.verifyEqual(E1, E0, 'beam centre');

end

function testVisualise(testCase)

  waist = 0.5;
  beam = ott.beam.GaussianDavis5(waist);
  
  h = figure();
  beam.visualise();
  axis image;
  close(h);
  
end
