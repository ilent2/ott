function tests = testDielectric
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../../');
end

function testConstructor(testCase)

  permittivity = 1.0;

  material = ott.beam.medium.Dielectric(permittivity);
  testCase.verifyEqual(material.relative_permittivity, permittivity, ...
      'permittivity');
  testCase.verifyEqual(material.relative_permeability, 1.0, 'permeability');
  
end

function testFromIndex(testCase)

  index = 1.0;

  material = ott.beam.medium.Dielectric.FromIndex(index);
  testCase.verifyEqual(material.index, index, ...
      'index');
  testCase.verifyEqual(material.relative_permeability, 1.0, 'permeability');

end
