function tests = testGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  beam = ott.beam.Gaussian();
  testCase.verifyEqual(beam.power, 1.0, 'power');
  testCase.verifyEqual(beam.mapping, 'sin', 'mapping');
  testCase.verifyEqual(beam.polfield, [1, 1i], 'polfield');
  testCase.verifyEqual(beam.polbasis, 'cartesian', 'polbasis');
  testCase.verifyEqual(beam.isGaussian, true, 'isGauss');

end

function testFromNa(testCase)

  NA = 0.9;
  beam = ott.beam.Gaussian.FromNa(NA);
  testCase.assertInstanceOf(beam, 'ott.beam.Gaussian');

end

function testIntensityMoment(testCase)

  waist = 15e-6;  % m
  power = 1.0;   % W
  beam = ott.beam.Gaussian('waist', waist, 'power', power, ...
    'polbasis', 'cartesian', 'polfield', [1;0]);
  beam.visNearfield('field', 'E2', 'range', [1,1]*6e-6)
  
%   bsc = beam.data;
%   bsc.power = power; % .* beam.impedance .* beam.wavenumber.^2;
  
%   E0 = beam.efieldRtp([0;0;0]).vxyz .* beam.wavenumber;
%   E0 = bsc.efieldRtp([0;0;0]./beam.wavelength).vxyz
%   H0 = bsc.hfieldRtp([0;0;0]).vxyz
  
  % Formula from Wikipedia
%   I0_trial = abs(E0(1)).^2 ./ (2*beam.impedance)
%   I0_target = 2*power ./ (pi * waist^2)
%   I0_trial / I0_target
  
%   testCase.verifyEqual(ints, beam1.power, 'RelTol', 1e-2, 'power');
%   testCase.verifyEqual(moment./beam1.speed, ...
%     beam1.force(emptyBeam), 'RelTol', 1e-2, 'AbsTol', 1e-23, 'moment');

end

