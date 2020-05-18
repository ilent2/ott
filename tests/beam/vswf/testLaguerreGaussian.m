function tests = testLaguerreGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructDefault(testCase)

  beam = ott.beam.vswf.LaguerreGaussian();

  testCase.verifyEqual(beam.waist, 1.0, 'waist');
  testCase.verifyEqual(beam.pmode, 0, 'azimuthal_mode');
  testCase.verifyEqual(beam.lmode, 0, 'radial_mode');

end

function testConstructOptional(testCase)

  waist = 0.1;
  azimuthal_mode = 1;
  radial_mode = 2;
  beam = ott.beam.vswf.LaguerreGaussian(waist, azimuthal_mode, radial_mode);

  testCase.verifyEqual(beam.waist, waist, 'waist');
  testCase.verifyEqual(beam.lmode, azimuthal_mode, 'azimuthal_mode');
  testCase.verifyEqual(beam.pmode, radial_mode, 'radial_mode');

end

function testModeErrors(testCase)

  import ott.beam.vswf.*;
  waist = 1.0;

  % Negative radial mode
  testCase.verifyError(@() LaguerreGaussian(waist, 0, -2), ...
    'ott:vswf:LaguerreGaussian:invalid_radial_mode');

  % non-integer radial mode
  testCase.verifyError(@() LaguerreGaussian(waist, 0, 1.5), ...
    'ott:vswf:LaguerreGaussian:invalid_radial_mode');

  % non-integer azimuthal mode
  testCase.verifyError(@() LaguerreGaussian(waist, 0.5, 0), ...
    'ott:vswf:LaguerreGaussian:invalid_azimuthal_mode');

end

