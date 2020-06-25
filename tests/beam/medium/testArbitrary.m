function tests = testArbitrary
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../../');
end

function testConstructor(testCase)

  permittivity = 1.0;
  permeability = 3.0;

  material = ott.beam.medium.Arbitrary(permittivity, permeability);
  testCase.verifyEqual(material.relative_permittivity, permittivity, ...
      'permittivity');
  testCase.verifyEqual(material.relative_permeability, permeability, ...
      'permeability');
  
end

