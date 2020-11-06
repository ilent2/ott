function tests = testEmpty
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  beam = ott.beam.Empty();
  testCase.verifyEqual(beam.position, [0;0;0], 'pos');
  testCase.verifyEqual(beam.rotation, eye(3), 'rot');
  testCase.verifyEqual(beam.index_medium, 1, 'index');
  testCase.verifyEqual(beam.omega, 2*pi*3e8/1064e-9, ...
      'RelTol', 1e-15, 'wavelength');

end

function testField(testCase)

  beam = ott.beam.Empty();
  xyz = randn(3, 5);
  target = zeros(size(xyz));
  
  testCase.verifyEqual(beam.efield(xyz), target, 'efield');
  testCase.verifyEqual(beam.hfield(xyz), target, 'hfield');
  testCase.verifyEqual(beam.ehfield(xyz), target, 'ehfield');
  
  testCase.verifyEqual(beam.efarfield(xyz), target, 'efarfield');
  testCase.verifyEqual(beam.hfarfield(xyz), target, 'hfarfield');
  testCase.verifyEqual(beam.ehfarfield(xyz), target, 'efarhfield');
  
  testCase.verifyEqual(beam.efieldRtp(xyz), target, 'efieldRtp');
  testCase.verifyEqual(beam.hfieldRtp(xyz), target, 'hfieldRtp');
  testCase.verifyEqual(beam.ehfieldRtp(xyz), target, 'ehfieldRtp');
  
  testCase.verifyEqual(beam.eparaxial(xyz), target, 'ep');
  testCase.verifyEqual(beam.hparaxial(xyz), target, 'hp');
  testCase.verifyEqual(beam.ehparaxial(xyz), target, 'ehp');
end
