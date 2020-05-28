function tests = testLaguerreGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  % Position arguments
  waist = 1.0;
  lmode = 1;
  pmode = 2;
  beam = ott.beam.abstract.LaguerreGaussian(waist, lmode, pmode);
  testCase.verifyEqual(beam.waist, waist);
  testCase.verifyEqual(beam.lmode, lmode);
  testCase.verifyEqual(beam.pmode, pmode);
  testCase.verifyEqual(beam.power, 1.0);

  % Named arguments
  beam = ott.beam.abstract.LaguerreGaussian('waist', waist, ...
      'lmode', lmode, 'pmode', pmode);
  testCase.verifyEqual(beam.waist, waist);
  testCase.verifyEqual(beam.lmode, lmode);
  testCase.verifyEqual(beam.pmode, pmode);
  testCase.verifyEqual(beam.power, 1.0);
end

function testConvertBeam(testCase)

  abs_beam = ott.beam.abstract.LaguerreGaussian(...
      'waist', 1.0, 'lmode', 5, 'pmode', 3);

  % VSWF Casts

  beam = ott.beam.vswf.Bsc(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.LaguerreGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.LaguerreGaussian, ...
      beam, abs_beam);

  beam = ott.beam.vswf.LaguerreGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.vswf.LaguerreGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.LaguerreGaussian, ...
      beam, abs_beam);

  % Paraxial

  beam = ott.beam.paraxial.Paraxial(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.LaguerreGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.LaguerreGaussian, ...
      beam, abs_beam);

  beam = ott.beam.paraxial.LaguerreGaussian(abs_beam);
  testCase.verifyClass(beam, 'ott.beam.paraxial.LaguerreGaussian');
  verifyProperties(testCase, ?ott.beam.abstract.LaguerreGaussian, ...
      beam, abs_beam);

end

