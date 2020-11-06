function tests = testBscBeam
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  index = 1.33;
  omega = 2*pi*3e8 ./ 532e-9;
  bsc = ott.bsc.Bsc([1;2;3], [4;5;6]);
  beam = ott.beam.BscBeam(bsc, 'index_medium', index, 'omega', omega);
  
  testCase.verifyEqual(beam.data, bsc, 'bsc');
  testCase.verifyEqual(beam.index_medium, index, 'index');
  testCase.verifyEqual(beam.omega, omega, 'omega');
end

function testGetData(testCase)

  bsc = ott.bsc.Bsc([1;2;3;4;5;6;7;8], zeros(8, 1));
  beam = ott.beam.BscBeam(bsc);
  
  % Original beam
  trial = ott.bsc.Bsc(beam);
  testCase.verifyEqual(trial, bsc, 'original');
  
  % Grow Nmax
  testCase.verifyError(@() ott.bsc.Bsc(beam, 10), ...
    'ott:beam:BscBeam:recalculate_not_implemented');
  
  % Rotate
  rbsc = bsc.rotateY(pi/2);
  rbeam = beam.rotateY(pi/2);
  trial = ott.bsc.Bsc(rbeam);
  testCase.verifyEqual(trial, rbsc, 'rotate');
  
  % Translate small
  tbsc = bsc.translateXyz([0.001;0;0], 'Nmax', 1);
  tbeam = beam.translateXyz([0.001;0;0]*beam.wavelength);
  trial = ott.bsc.Bsc(tbeam);
  testCase.verifyEqual(trial, tbsc, 'translate');

end

function testNearfieldFunctions(testCase)

  index = 1.33;
  omega = 2*pi*3e8 ./ 532e-9;
  bsc = ott.bsc.Bsc([1;2;3], [4;5;6]);
  beam = ott.beam.BscBeam(bsc, 'index_medium', index, 'omega', omega);
  k = beam.wavenumber;
  kZ = beam.wavenumber * beam.impedance;
  
  xyz = randn(3, 5);
  rtp = ott.utils.xyz2rtp(xyz);
  
  targetE = bsc.efieldRtp(rtp./[beam.wavelength;1;1])./k;
  targetH = bsc.hfieldRtp(rtp./[beam.wavelength;1;1])./kZ;
  
  trialE = beam.efieldRtp(rtp);
  testCase.verifyEqual(trialE.vrtp, targetE.vrtp, 'RelTol', 1e-15, 'efieldRtp');
  trialE = beam.efield(xyz);
  testCase.verifyEqual(trialE.vrtp, targetE.vrtp, 'RelTol', 1e-15, 'efield');
  
  trialH = beam.hfieldRtp(rtp);
  testCase.verifyEqual(trialH.vrtp, targetH.vrtp, 'RelTol', 1e-15, 'hfieldRtp');
  trialH = beam.hfield(xyz);
  testCase.verifyEqual(trialH.vrtp, targetH.vrtp, 'RelTol', 1e-15, 'hfield');
  
  [trialE, trialH] = beam.ehfieldRtp(rtp);
  testCase.verifyEqual(trialE.vrtp, targetE.vrtp, 'RelTol', 1e-15, 'ehfieldRtp e');
  testCase.verifyEqual(trialH.vrtp, targetH.vrtp, 'RelTol', 1e-15, 'ehfieldRtp h');
  [trialE, trialH] = beam.ehfield(xyz);
  testCase.verifyEqual(trialE.vrtp, targetE.vrtp, 'RelTol', 1e-15, 'ehfield e');
  testCase.verifyEqual(trialH.vrtp, targetH.vrtp, 'RelTol', 1e-15, 'ehfield h');
  
end

function testFarfieldFunctions(testCase)

  index = 1.33;
  omega = 2*pi*3e8 ./ 532e-9;
  bsc = ott.bsc.Bsc([1;2;3], [4;5;6]);
  beam = ott.beam.BscBeam(bsc, 'index_medium', index, 'omega', omega);
  k = beam.wavenumber;
  Z = beam.impedance;
  
  rtp = randn(3, 5);

  targetE = bsc.efarfield(rtp)./k;
  targetH = ott.utils.FieldVectorSph(...
          -1i .* targetE.vrtp([1, 3, 2], :)./Z, rtp);
  
  trialE = beam.efarfield(rtp);
  testCase.verifyEqual(trialE.vrtp, targetE.vrtp, 'RelTol', 1e-15, 'e');
  
  trialH = beam.hfarfield(rtp);
  testCase.verifyEqual(trialH.vrtp, targetH.vrtp, 'RelTol', 1e-15, 'h');
  
  [trialE, trialH] = beam.ehfarfield(rtp);
  testCase.verifyEqual(trialE.vrtp, targetE.vrtp, 'RelTol', 1e-15, 'eh e');
  testCase.verifyEqual(trialH.vrtp, targetH.vrtp, 'RelTol', 1e-15, 'eh h');
  
end

function testForceFunctions(testCase)

  bsc1 = ott.bsc.Bsc([1;2;3], 0);
  bsc2 = ott.bsc.Bsc([3;2;1], 0);
  
  beam1 = ott.beam.BscBeam(bsc1);
  beam2 = ott.beam.BscBeam(bsc2);
  
  target = bsc1.force(bsc2) ./ beam1.speed;
  trial = beam1.force(beam2);
  testCase.verifyEqual(target, trial, 'force');
  
  target = bsc1.torque(bsc2) ./ beam1.omega;
  trial = beam1.torque(beam2);
  testCase.verifyEqual(target, trial, 'torque');
  
  target = bsc1.spin(bsc2) ./ beam1.omega;
  trial = beam1.spin(beam2);
  testCase.verifyEqual(target, trial, 'spin');

end
