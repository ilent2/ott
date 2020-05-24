function tests = testGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  % Position arguments
  waist = 0.7;
  beam = ott.beam.abstract.Gaussian(waist);
  testCase.verifyEqual(beam.waist, waist);
  testCase.verifyEqual(beam.power, 1.0);

  % Named arguments
  waist = 0.3;
  beam = ott.beam.abstract.Gaussian('waist', waist);
  testCase.verifyEqual(beam.waist, waist);
  testCase.verifyEqual(beam.power, 1.0);
end

function testConvertBeam(testCase)

  waist = 1.0;
  abs_beam = ott.beam.abstract.Gaussian(waist);

  % Beam casts

  beam = ott.beam.Beam(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.GaussianDavis5');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.GaussianDavis5(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.GaussianDavis5');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  % VSWF Casts

  beam = ott.beam.vswf.Bsc(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.Gaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

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

  beam = ott.beam.paraxial.HermiteGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.HermiteGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.paraxial.LaguerreGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.LaguerreGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.paraxial.InceGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.InceGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);
end

