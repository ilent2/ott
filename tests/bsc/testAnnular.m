function tests = testAnnular
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testFromBessel(testCase)

  Nmax = 3;
  theta = 0.5;
  Etp = [1;0];
  lmode = 0;
  beam = ott.bsc.Annular.FromBessel(Nmax, theta, Etp, lmode);

  testCase.assertClass(beam, 'ott.bsc.Annular');
  testCase.verifyEqual(beam.Nmax, Nmax, 'Nmax');
  testCase.verifyEqual(beam.theta, theta, 'theta');

end

function testFromMathieu(testCase)

  theta = pi/4;
  morder = 1;
  ellipticity = 1;
  parity = 'even';
  Nmax = 5;
  [beam, besselWeights] = ott.bsc.Annular.FromMathieu(...
      theta, morder, ellipticity, parity, ...
      'Nmax', Nmax);

  testCase.assertClass(beam, 'ott.bsc.Annular');
  testCase.verifyEqual(beam(1).Nmax, Nmax, 'Nmax');
  testCase.verifyEqual(beam(1).theta, theta, 'theta');
  testCase.verifySize(beam, [1, numel(besselWeights)], 'size');

end

function testFromWebber(testCase)

  theta = pi/4;
  alpha = 1;
  parity = 'even';
  Nmax = 5;
  [beam, besselWeights] = ott.bsc.Annular.FromWebber(...
      theta, alpha, parity, 'Nmax', Nmax);

  testCase.assertClass(beam, 'ott.bsc.Annular');
  testCase.verifyEqual(beam(1).Nmax, Nmax, 'Nmax');
  testCase.verifyEqual(beam(1).theta, theta, 'theta');
  testCase.verifySize(beam, [1, numel(besselWeights)], 'size');

end

function testArray(testCase)

  Nmax = 3;
  theta = [0.5, 0.2];
  Etp = [1;0];
  lmode = 0;
  beam = ott.bsc.Annular.FromBessel(Nmax, theta, Etp, lmode);

  testCase.assertClass(beam, 'ott.bsc.Annular');
  testCase.assertSize(beam, [1, 2], 'size');
  testCase.verifyEqual([beam.Nmax], [1,1]*Nmax, 'Nmax');
  testCase.verifyEqual([beam.theta], theta, 'theta');
end

function testTranslateZ(testCase)
  % Test translation in comparision to normal bsc translation

  Nmax = 40;
  theta = 0.5;
  Etp = [1;0];
  lmode = 1;
  beam = ott.bsc.Annular.FromBessel(Nmax, theta, Etp, lmode);

  % Calculate translation with Annular
  tbeam = beam.translateZ(1);

  % Calculate translation with Bsc
  bsc = ott.bsc.Bsc(beam);
  tbsc = bsc.translateZ(1);

  % Compare coefficients
  ab1 = tbeam.getCoefficients(1:ott.utils.combined_index(10, 10));
  ab2 = tbsc.getCoefficients(1:ott.utils.combined_index(10, 10));
  testCase.verifyEqual(ab1, ab2, 'RelTol', 1.0e-6);

end
