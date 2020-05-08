function tests = testVacuum
  % Test for beam combination functionality
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testWrongArgCount(testCase)

  testCase.verifyError(@() ott.beam.medium.Vacuum(), ...
    'ott:beam:medium:Vacuum:wrong_arg_count');
  
end

function testSi(testCase)

  vacuum = ott.beam.medium.Vacuum.SiUnits;
  testCase.verifyEqual(vacuum.permittivity, 8.854e-12, 'RelTol', 1e-3);
  testCase.verifyEqual(vacuum.permeability, 1.257e-6, 'RelTol', 1e-3);
  testCase.verifyEqual(vacuum.speed, 2.998e8, 'RelTol', 1e-3);
  
end

function testUnit(testCase)

  vacuum = ott.beam.medium.Vacuum.Unitary;
  testCase.verifyEqual(vacuum.permittivity, 1.0, 'RelTol', 1e-3);
  testCase.verifyEqual(vacuum.permeability, 1.0, 'RelTol', 1e-3);
  testCase.verifyEqual(vacuum.speed, 1.0, 'RelTol', 1e-3);
  
end
