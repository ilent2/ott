function tests = testBsc
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

  testCase.verifyClass(beam, 'ott.bsc.Annular');
  testCase.verifyEqual(beam.Nmax, Nmax, 'Nmax');
  testCase.verifyEqual(beam.theta, theta, 'Nmax');

end

function testTranslateZ(testCase)
  % Test translation in comparision to normal bsc translation

  Nmax = 40;
  theta = 0.5;
  Etp = [1;0];
  lmode = 1;
  beam = ott.bsc.Annular.FromBessel(Nmax, theta, Etp, lmode);

  % Calculate translation with Annular
  tbeam = beam.translateZInternal(1);

  % Calculate translation with Bsc
  bsc = ott.bsc.Bsc(beam);
  tbsc = bsc.translateZInternal(1);

  % Compare coefficients
  ab1 = tbeam.getCoefficients(1:ott.utils.combined_index(10, 10));
  ab2 = tbsc.getCoefficients(1:ott.utils.combined_index(10, 10));
  testCase.verifyEqual(ab1, ab2, 'RelTol', 1.0e-6);

end
