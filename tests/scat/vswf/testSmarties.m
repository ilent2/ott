function tests = testTmatrixSmarties
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../');
end

function testMie(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  Tmie = ott.TmatrixMie(0.1, 'index_relative', 1.2);

  Ttest = ott.TmatrixSmarties.simple('sphere', 0.1, ...
    'index_relative', 1.2);

  testCase.verifyThat(Ttest.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix');

  tol = 1.0e-2;
  Tebcm_data = Ttest.data ./ max(abs(Ttest.data(:)));
  Tmie_data = Tmie.data ./ max(abs(Tmie.data(:)));
  testCase.verifyThat(Tebcm_data, IsEqualTo(Tmie_data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance');

end

function testMie2(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  Tmie = ott.TmatrixMie(0.1e-6, 'index_medium', 1.0, ...
    'index_particle', 1.2, 'wavelength0', 1e-6);

  Ttest = ott.TmatrixSmarties.simple('sphere', 0.1e-6, ...
    'index_medium', 1.0, 'index_particle', 1.2, 'wavelength0', 1e-6);

  testCase.verifyThat(Ttest.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix');

  tol = 1.0e-2;
  Tebcm_data = Ttest.data ./ max(abs(Ttest.data(:)));
  Tmie_data = Tmie.data ./ max(abs(Tmie.data(:)));
  testCase.verifyThat(Tebcm_data, IsEqualTo(Tmie_data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance');

end

function testEllipse(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  Tebcm = ott.TmatrixEbcm.simple('ellipsoid', [0.1, 0.1, 0.09], ...
    'index_relative', 1.2);

  Ttest = ott.TmatrixSmarties.simple('ellipsoid', [0.1, 0.1, 0.09], ...
    'index_relative', 1.2);

  testCase.verifyThat(Ttest.Nmax, IsEqualTo(Tebcm.Nmax), ...
      'Nmax does not match EBCM T-matrix');

  tol = 1.0e-2;
  Tebcm_data = Ttest.data ./ max(abs(Ttest.data(:)));
  Tmie_data = Tebcm.data ./ max(abs(Tebcm.data(:)));
  testCase.verifyThat(Tebcm_data, IsEqualTo(Tmie_data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match EBCM T-matrix within tolerance');

end

function testInternal(testCase)
  % Test construction of internal T-matrix
  
  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  Tmie = ott.TmatrixMie(0.1, 'index_relative', 1.2, 'internal', true);

  Ttest = ott.TmatrixSmarties.simple('sphere', 0.1, ...
      'index_relative', 1.2, 'internal', true);

  testCase.verifyThat(Ttest.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix');

  tol = 1.0e-2;
  Tebcm_data = Ttest.data ./ max(abs(Ttest.data(:)));
  Tmie_data = Tmie.data ./ max(abs(Tmie.data(:)));
  testCase.verifyThat(Tebcm_data, IsEqualTo(Tmie_data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance');
    
end