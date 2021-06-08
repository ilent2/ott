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
    
  % Set methods
  beam.speed = 2;
  testCase.verifyEqual(beam.index_medium, 3e8/2, 'change speed');
  testCase.verifyError(@setWavenumber, 'ott:beam:Beam:set_wavenumber', 'k');
  testCase.verifyError(@setWavelength, 'ott:beam:Beam:set_wavelength', 'w');
  a = beam.setWavenumber(1, 'fixedSpeed');
  a = beam.setWavelength(2, 'fixedFrequency');
  
  function setWavenumber()
    beam.wavenumber = 1;
  end

  function setWavelength()
    beam.wavelength = 1;
  end

end

function testForce(testCase)

  b1 = ott.beam.Empty();
  b2 = ott.beam.Empty();
  
  f = b1.force(b2);
  testCase.verifyEqual(f, [0;0;0], 'force');
  
  f = b1.torque(b2);
  testCase.verifyEqual(f, [0;0;0], 'torque');
  
  f = b1.spin(b2);
  testCase.verifyEqual(f, [0;0;0], 'spin');

end

function testField(testCase)

  beam = ott.beam.Empty();
  xyz = randn(3, 5);
  target = ott.utils.FieldVectorCart(zeros(size(xyz)), xyz);
    
  rtp = mod(randn(3, 5), pi);
  targetS = ott.utils.FieldVectorSph(zeros(size(rtp)), rtp);
  targetF = ott.utils.FieldVectorSph(zeros(size(rtp)), ...
      [ones(1, size(rtp, 2)); rtp(2:3, :)]);
    
  xy = randn(2, 5);
  
  testCase.verifyEqual(beam.efield(xyz), target, 'efield');
  testCase.verifyEqual(beam.hfield(xyz), target, 'hfield');
  testCase.verifyEqual(beam.ehfield(xyz), target, 'ehfield');
  
  testCase.verifyEqual(beam.efarfield(rtp), targetF, 'efarfield');
  testCase.verifyEqual(beam.hfarfield(rtp), targetF, 'hfarfield');
  testCase.verifyEqual(beam.ehfarfield(rtp), targetF, 'efarhfield');
  
  testCase.verifyEqual(beam.efieldRtp(rtp), targetS, 'efieldRtp');
  testCase.verifyEqual(beam.hfieldRtp(rtp), targetS, 'hfieldRtp');
  testCase.verifyEqual(beam.ehfieldRtp(rtp), targetS, 'ehfieldRtp');
  
  testCase.verifyEqual(beam.eparaxial(xy).vrtp(), zeros(size(xyz)), 'ep');
  testCase.verifyEqual(beam.hparaxial(xy).vrtp(), zeros(size(xyz)), 'hp');
  testCase.verifyEqual(beam.ehparaxial(xy).vrtp(), zeros(size(xyz)), 'ehp');
  
  testCase.verifyEqual(beam.intensityMoment(), zeros(3, 1), 'im');
end
