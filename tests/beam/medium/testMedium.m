function tests = testMedium
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../../');
end

function testConstructor(testCase)

  water = ott.beam.medium.Generic.Water;
  vacuum = ott.beam.medium.Vacuum.BaseSi;
  freq = 1.0;
  
  medium = ott.beam.medium.Medium(water, freq, vacuum);
  
  testCase.assertClass(medium, 'ott.beam.medium.Medium', 'mat');
  
  % Test properties
  testCase.verifyEqual(medium.material, water, 'mat');
  testCase.verifyEqual(medium.vacuum, vacuum, 'vacuum');
  testCase.verifyEqual(medium.frequency, freq, 'freq');
  
  % Test dependent properties
  testCase.verifyEqual(medium.permittivity, ...
    water.relative_permittivity .* vacuum.permittivity, 'permittivity');
  testCase.verifyEqual(medium.permeability, ...
    water.relative_permeability .* vacuum.permeability, 'permeability');
  testCase.verifyEqual(medium.speed, ...
    vacuum.speed ./ water.index, 'speed');
  testCase.verifyEqual(medium.index, ...
    water.index, 'index');
  testCase.verifyEqual(medium.impedance, ...
    sqrt(medium.permeability ./ medium.permittivity), 'impedance');
  
end

function testDefaultFrequency(testCase)

  % Clear persistent variables and test default
  clear ott.beam.medium.Medium;
  testCase.verifyEqual(ott.beam.medium.Medium.DefaultFrequency, 2*pi, 'default');
  
  % Test assign
  val = 1;
  ott.beam.medium.Medium.DefaultFrequency(val);
  testCase.verifyEqual(ott.beam.medium.Medium.DefaultFrequency, val, 'getset');
end

function testDefaultVacuum(testCase)

  % Clear persistent variables and test default
  clear ott.beam.medium.Medium;
  testCase.verifyEqual(ott.beam.medium.Medium.DefaultVacuum, ...
    ott.beam.medium.Vacuum.Unitary, 'default');
  
  % Test assign
  val = ott.beam.medium.Vacuum.BaseSi;
  ott.beam.medium.Medium.DefaultVacuum(val);
  testCase.verifyEqual(ott.beam.medium.Medium.DefaultVacuum, val, 'getset');
end
