function tests = testGaussian
  % Test for Gaussian and the base Paraxial class casts
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

function testGaussianConvertBeam(testCase)

  waist = 1.0;
  abs_beam = ott.beam.abstract.Gaussian(waist);

  beam = ott.beam.vswf.Bsc(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.Gaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);

  beam = ott.beam.paraxial.HermiteGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.HermiteGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);
  
end

