function tests = testGaussian
  % Test for Gaussian and the base Paraxial class casts
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testParaxialConvertBeam(testCase)

  waist = 1.0;
  abs_beam = ott.beam.abstract.Gaussian(waist);

  % Beam casts

  beam = ott.beam.Beam(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.Gaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.GaussianDavis5(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.GaussianDavis5');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  % VSWF Casts

  beam = ott.beam.vswf.Gaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.Gaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.vswf.HermiteGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.HermiteGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.vswf.LaguerreGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.LaguerreGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.vswf.InceGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.InceGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  % Paraxial casts

  beam = ott.beam.paraxial.Gaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.Gaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.paraxial.LaguerreGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.LaguerreGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.paraxial.HermiteGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.HermiteGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  % Abstract casts

  beam = ott.beam.abstract.Gaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.abstract.Gaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.abstract.LaguerreGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.abstract.LaguerreGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.abstract.HermiteGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.abstract.HermiteGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.abstract.InceGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.abstract.InceGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);
end
