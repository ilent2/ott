function tests = testPolarizabilityLdr
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testArgumentCount(testCase)

  spacing = 1.0;
  index = 1.0;
  kvec = [0, 0, 1];
  E0 = [1, 0, 0];

  alpha = ott.utils.polarizability.LDR(spacing,index,kvec,E0);
  
  alpha = ott.utils.polarizability.LDR(spacing,index,kvec,E0, 'k0', 2*pi);
  
  alpha = ott.utils.polarizability.LDR(spacing,index);
  
  alpha = ott.utils.polarizability.LDR(spacing,index, 'k0', 2*pi);

end

function testAlphaValues(testCase)

  spacing = 1.0;
  index = [1, 2];
  
  alpha = ott.utils.polarizability.LDR(spacing,index);
  
  targetAlpha = [0; -0.002627924391542 + 0.004518926661470i];
  
  testCase.verifyEqual(alpha, targetAlpha, ...
    'AbsTol', 1e-15, ...
    'Alpha doesnt match target values');
end
