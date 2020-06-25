function tests = testRelative
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../../');
end

function testConstructor(testCase)

  m1 = ott.beam.medium.Generic.Polystyrene;
  m2 = ott.beam.medium.Generic.Water;
  
  % Constructor is called by operator
  mat = m1 ./ m2;
  testCase.verifyEqual(mat.material1, m1, 'm1');
  testCase.verifyEqual(mat.material2, m2, 'm2');
  
end
