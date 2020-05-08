function tests = testUnmatchedArgs
  % Test Cartesian and Spherical coordinate transform functions
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testMultipleArgs(testCase)

  p = inputParser;
  p.KeepUnmatched = true;
  p.addParameter('test', []);
  p.addParameter('more', []);
  p.parse('hah', 1, 'nah', 4, 'lay', 3);
  
  unmatched = ott.utils.unmatchedArgs(p);
  
  testCase.verifyEqual(size(unmatched), [2, 3], 'sz');
  
  p = inputParser;
  p.addParameter('hah', []);
  p.addParameter('nah', []);
  p.addParameter('lay', []);
  p.parse(unmatched{:});
  
  unmatched = ott.utils.unmatchedArgs(p);
  testCase.verifyEqual(numel(unmatched), 0, 'empty');
  
end
