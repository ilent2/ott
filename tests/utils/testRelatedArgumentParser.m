function tests = testRelatedArgumentParser
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testExample(testCase)

  p = ott.utils.RelatedArgumentParser;
  p.addRequired('index_medium', 1.0);
  p.addOptional('speed0', 2.0);
  p.addOptional('speed');
  p.addRule('index_medium', @(s1, s0) s0 ./ s1, 'speed', 'speed0');
  p.parse('speed', 1.0);
  
  testCase.verifyEqual(p.RequiredResults.index_medium, 2, ...
    'AbsTol', 1e-15);

end

function testAddRuleErrors(testCase)

  p = ott.utils.RelatedArgumentParser;
  p.addRequired('index_medium', 1.0);
  p.addOptional('speed0', 2.0);
  p.addOptional('speed');
  
  testCase.verifyError(...
    @() p.addRule('not_a_var', @() [], 'speed'), ...
    'ott:utils:RelatedArgumentParser:var_not_a_required');
  
  testCase.verifyError(...
    @() p.addRule('index_medium', @() [], 'not_a_param'), ...
    'ott:utils:RelatedArgumentParser:arg_not_a_parameter');

end

