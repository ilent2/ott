function tests = testInceGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  % Position arguments
  waist = 1.0;
  lmode = 3;
  porder = 5;
  parity = 'even';
  ellip = 1.0;
  beam = ott.beam.abstract.InceGaussian(waist, lmode, porder, parity, ellip);
  testCase.verifyEqual(beam.waist, waist);
  testCase.verifyEqual(beam.lmode, lmode);
  testCase.verifyEqual(beam.porder, porder);
  testCase.verifyEqual(beam.parity, parity);
  testCase.verifyEqual(beam.ellipticity, ellip);
  testCase.verifyEqual(beam.power, 1.0);

  % Named arguments
  beam = ott.beam.abstract.InceGaussian('waist', waist, 'lmode', lmode, ...
      'porder', porder, 'parity', parity, 'ellipticity', ellip);
  testCase.verifyEqual(beam.waist, waist);
  testCase.verifyEqual(beam.lmode, lmode);
  testCase.verifyEqual(beam.porder, porder);
  testCase.verifyEqual(beam.parity, parity);
  testCase.verifyEqual(beam.ellipticity, ellip);
  testCase.verifyEqual(beam.power, 1.0);
end

function testConvertBeam(testCase)

  beam = ott.beam.abstract.InceGaussian('waist', 1.0, 'lmode', 1, ...
      'porder', 5, 'parity', 'even', 'ellipticity', 1.0);

  % Beam casts

  beam = ott.beam.Beam(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.InceGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.InceGaussian, ...
      beam, abs_beam);

  % VSWF Casts

  beam = ott.beam.vswf.Bsc(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.InceGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.InceGaussian, ...
      beam, abs_beam);

  beam = ott.beam.vswf.InceGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.InceGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.InceGaussian, ...
      beam, abs_beam);
end

function testGaussianCasts(testCase)

  beam = ott.beam.abstract.InceGaussian('waist', 1.0, 'lmode', 0, ...
      'porder', 0, 'parity', 'even', 'ellipticity', 1.0);

  % VSWF casts

  beam = ott.beam.vswf.Gaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.Gaussian');
  verifyProperties(testCase, ?ott.beam.abstract.InceGaussian, ...
      beam, abs_beam);

  beam = ott.beam.vswf.LaguerreGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.LaguerreGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.InceGaussian, ...
      beam, abs_beam);

  beam = ott.beam.vswf.HermiteGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.HermiteGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.InceGaussian, ...
      beam, abs_beam);

  % Paraxial casts

  beam = ott.beam.paraxial.LaguerreGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.LaguerreGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.InceGaussian, ...
      beam, abs_beam);

  beam = ott.beam.paraxial.HermiteGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.HermiteGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.InceGaussian, ...
      beam, abs_beam);
end

