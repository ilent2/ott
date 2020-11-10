function tests = testStrata
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  depths = [1;2];
  normal = [0;0;1];
  shape = ott.shape.Strata(depths, normal);

  testCase.verifyEqual(shape.normal, normal, 'normal');
  testCase.verifyEqual(shape.depths, depths, 'depths');
  
  % Scale
  shape = shape * 2;
  testCase.verifyEqual(shape.depths, 2*depths, 'S depths');
  testCase.verifyEqual(shape.normal, normal, 'S normals');
  

end

