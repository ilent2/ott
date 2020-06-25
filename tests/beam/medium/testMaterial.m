function tests = testMaterial
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../../');
end

% Abstract class, no constructor

function testProperties(testCase)

  material = ott.beam.medium.Generic.Water;
  
  testCase.verifyEqual(material.index, ...
    sqrt(material.relative_permittivity .* material.relative_permeability), 'index');
  testCase.verifyEqual(material.isIsotropic, true, 'isotropic');
  testCase.verifyEqual(material.isConductive, false, 'conductive');
  
end
