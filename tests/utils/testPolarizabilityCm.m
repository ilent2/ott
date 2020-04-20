function tests = testPolarizabilityCm
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testAlphaValues(testCase)

  spacing = 1.0;
  index = [1, 2];
  
  alpha = ott.utils.polarizability.CM(spacing,index);
  
  targetAlpha = [0; 0.119366207318922];
  
  testCase.verifyEqual(alpha, targetAlpha, ...
    'AbsTol', 1e-15, ...
    'Alpha doesnt match target values');
end
