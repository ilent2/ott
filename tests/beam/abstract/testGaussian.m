function tests = testGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  waist = 1.0;
  beam = ott.beam.abstract.Gaussian(waist);
  testCase.verifyEqual(beam.waist, waist);
end

function verifyProperties(testCase, meta, newbeam, oldbeam)
  props = {meta.PropertyList(~[meta.PropertyList.Dependent]).Name};
  for p = props
    testCase.verifyEqual(newbeam.(p{1}), oldbeam.(p{1}), p{1});
  end
end

function testConvertBeam(testCase)

  waist = 1.0;
  abs_beam = ott.beam.abstract.Gaussian(waist);
  beam = ott.beam.Beam(abs_beam);
  
  testCase.verifyClass(beam, 'ott.beam.GaussianDavis5');
  verifyProperties(testCase, ?ott.beam.abstract.Gaussian, beam, abs_beam);
  
end

