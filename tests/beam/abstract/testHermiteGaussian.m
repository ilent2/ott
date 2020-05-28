function tests = testHermiteGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  % Position arguments
  waist = 1.0;
  mmode = 1;
  nmode = 2;
  beam = ott.beam.abstract.HermiteGaussian(waist, mmode, nmode);
  testCase.verifyEqual(beam.waist, waist);
  testCase.verifyEqual(beam.mmode, mmode);
  testCase.verifyEqual(beam.nmode, nmode);
  testCase.verifyEqual(beam.power, 1.0);

  % Named arguments
  beam = ott.beam.abstract.HermiteGaussian('waist', waist, ...
      'mmode', mmode, 'nmode', nmode);
  testCase.verifyEqual(beam.waist, waist);
  testCase.verifyEqual(beam.mmode, mmode);
  testCase.verifyEqual(beam.nmode, nmode);
  testCase.verifyEqual(beam.power, 1.0);
end

function testConvertBeam(testCase)

  abs_beam = ott.beam.abstract.HermiteGaussian(...
      'waist', 1.0, 'mmode', 5, 'nmode', 3);

  % VSWF Casts

  beam = ott.beam.vswf.Bsc(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.HermiteGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.HermiteGaussian, ...
      beam, abs_beam);

  beam = ott.beam.vswf.HermiteGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.HermiteGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.HermiteGaussian, ...
      beam, abs_beam);

  % Paraxial

  beam = ott.beam.paraxial.Paraxial(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.HermiteGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.HermiteGaussian, ...
      beam, abs_beam);

  beam = ott.beam.paraxial.HermiteGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.HermiteGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.HermiteGaussian, ...
      beam, abs_beam);

end


