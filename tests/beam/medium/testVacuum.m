function tests = testVacuum
  % Test for beam combination functionality
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruct(testCase)

  permittivity = 8.854e-12;
  permeability = 1.257e-6;
  
  vacuum = ott.beam.medium.Vacuum(permittivity, permeability);
  
  testCase.verifyEqual(vacuum.permittivity, permittivity, 'permittivity');
  testCase.verifyEqual(vacuum.permeability, permeability, 'permeability');
  
  testCase.verifyEqual(vacuum.speed, 1./sqrt(permittivity*permeability), 'speed');
  
end

function testCastMaterial(testCase)

  vacuum = ott.beam.medium.Vacuum.Unitary;
  material = ott.beam.medium.Material(vacuum);
  
  testCase.verifyTrue(isa(material, 'ott.beam.medium.Material'), 'type');
end

function testCastMedium(testCase)

  vacuum = ott.beam.medium.Vacuum.Unitary;
  medium = ott.beam.medium.Medium(vacuum);
  
  testCase.verifyTrue(isa(medium, 'ott.beam.medium.Medium'), 'type');
end
