function tests = testScattered
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstructDefault(testCase)

  beam = ott.beam.Scattered();
  testCase.verifyEqual(beam.scattered, ott.beam.Empty(), 'scattered');
  testCase.verifyEqual(beam.incident, ott.beam.Empty(), 'incident');
  testCase.verifyEqual(beam.particle, [], 'particle');
  testCase.verifyEqual(beam.internal, [], 'internal');
  
  beam.defaultVisRangeInternal();

end

function testFieldRtp(testCase)

  beam = ott.beam.Scattered();
  rtp = [zeros(3, 1), [10;0;0]];
  xyz = ott.utils.rtp2xyz(rtp);
  
  target = zeros(3, 2);
  testCase.verifyEqual(beam.efieldRtp(rtp), target, 'er no shape');
  testCase.verifyEqual(beam.efield(xyz), target, 'e no shape');
  testCase.verifyEqual(beam.hfieldRtp(rtp), target, 'hr no shape');
  testCase.verifyEqual(beam.hfield(xyz), target, 'h no shape');
  testCase.verifyEqual(beam.ehfieldRtp(rtp), target, 'ehr no shape');
  testCase.verifyEqual(beam.ehfield(xyz), target, 'eh no shape');
  
  beam.shape = ott.shape.Sphere(1);
  target = [nan(3, 1), zeros(3, 1)];
  testCase.verifyEqual(beam.efieldRtp(rtp), target, 'er shape');
  testCase.verifyEqual(beam.efield(xyz), target, 'e shape');
  testCase.verifyEqual(beam.hfieldRtp(rtp), target, 'hr shape');
  testCase.verifyEqual(beam.hfield(xyz), target, 'h shape');
  testCase.verifyEqual(beam.ehfieldRtp(rtp), target, 'ehr shape');
  testCase.verifyEqual(beam.ehfield(xyz), target, 'eh shape');
  
end

function testFarfield(testCase)

  beam = ott.beam.Scattered();
  rtp = [10;0;0];
  xy = [0;0];

  % Test coverage
  beam.efarfield(rtp);
  beam.hfarfield(rtp);
  beam.ehfarfield(rtp);
  beam.hparaxial(xy);
  beam.eparaxial(xy);
  beam.ehparaxial(xy);

end

function testForceTorque(testCase)

  beam = ott.beam.Scattered();
  rtp = [10;0;0];
  xy = [0;0];
  
  beam.force();
  beam.torque();
  beam.spin();
  beam.forcetorque();
  
end
