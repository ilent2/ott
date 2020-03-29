function tests = testBscPmGauss
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../');
end

function testConstruct(testCase)

  beam = ott.BscPmGauss();

  import matlab.unittest.constraints.IsEqualTo;
  testCase.verifyThat(beam.Nbeams, IsEqualTo(1), ...
    'Incorrect number of beams stored');

  testCase.verifyThat(beam.mode, IsEqualTo([0 0]), ...
    'Incorrect number of beams stored');

end

function testLgModeErrors(testCase)

  % Negative radial mode
  testCase.verifyError(@() ott.BscPmGauss('lg', [-2, 0]), ...
    'ott:BscPmGauss:invalid_radial_mode');

  % non-integer radial mode
  testCase.verifyError(@() ott.BscPmGauss('lg', [1.5, 0]), ...
    'ott:BscPmGauss:invalid_radial_mode');

  % non-integer azimuthal mode
  testCase.verifyError(@() ott.BscPmGauss('lg', [0, 0.5]), ...
    'ott:BscPmGauss:invalid_azimuthal_mode');

  % Too many mode numbers
  testCase.verifyError(@() ott.BscPmGauss('lg', [0, 0, 5]), ...
    'ott:BscPmGauss:wrong_mode_length');

end

function testHgModeErrors(testCase)

  % Too many mode numbers
  testCase.verifyError(@() ott.BscPmGauss('hg', [0, 0, 5]), ...
    'ott:BscPmGauss:wrong_mode_length');

end

function testIgModeErrors(testCase)

  % Too many mode numbers
  testCase.verifyError(@() ott.BscPmGauss('ig', [0, 0, 5]), ...
    'ott:BscPmGauss:wrong_mode_length');

end

function testAngularScaling(testCase)

  beam = ott.BscPmGauss('angular_scaling', 'sintheta');
  testCase.verifyEqual(beam.angular_scaling, 'sintheta', ...
    'sintheta Angular scaling not set correctly');

  beam = ott.BscPmGauss('angular_scaling', 'tantheta');
  testCase.verifyEqual(beam.angular_scaling, 'tantheta', ...
    'tantheta Angular scaling not set correctly');
  
end
