function tests = testPolarizabilityLdr
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testAlphaValues(testCase)

  spacing = 1.0;
  index = [1, 2];
  
  alpha = ott.utils.polarizability.FCD(spacing,index);
  
  targetAlpha = [0; 0.000029131505380 - 0.003023300651752i];
  
  testCase.verifyEqual(alpha, targetAlpha, ...
    'AbsTol', 1e-15, ...
    'Alpha doesnt match target values');
end
