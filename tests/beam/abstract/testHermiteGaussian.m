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
  waist = 1.0;
  lmode = 1;
  pmode = 2;
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

  % Beam casts

  beam = ott.beam.Beam(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.HermiteGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.HermiteGaussian, ...
      beam, abs_beam);

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

function testGaussianCasts(testCase)

  abs_beam = ott.beam.abstract.HermiteGaussian(...
      'waist', 1.0, 'nmode', 0, 'mmode', 0);

  % VSWF casts

  beam = ott.beam.vswf.Gaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.Gaussian');
  verifyProperties(testCase, ?ott.beam.abstract.HermiteGaussian, ...
      beam, abs_beam);

  beam = ott.beam.vswf.LaguerreGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.LaguerreGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.HermiteGaussian, ...
      beam, abs_beam);

  beam = ott.beam.vswf.InceGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.InceGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.HermiteGaussian, ...
      beam, abs_beam);

  % Paraxial casts

  beam = ott.beam.paraxial.LaguerreGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.LaguerreGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.HermiteGaussian, ...
      beam, abs_beam);

  beam = ott.beam.paraxial.InceGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.InceGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.HermiteGaussian, ...
      beam, abs_beam);
end

