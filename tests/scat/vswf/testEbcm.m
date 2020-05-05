function tests = tmatrixebcm
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  Tebcm = ott.TmatrixEbcm.simple('sphere', 1.0, 'index_relative', 1.2);

  Tmie = ott.TmatrixMie(1.0, 'index_relative', 1.2);

  testCase.verifyThat(Tebcm.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix');

  tol = 1.0e-2;
  Tebcm_data = Tebcm.data ./ max(abs(Tebcm.data(:)));
  Tmie_data = Tmie.data ./ max(abs(Tmie.data(:)));
  testCase.verifyThat(Tebcm_data, IsEqualTo(Tmie_data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance');

end

function testInternal(testCase)
  % Test construction of internal T-matrix
  
  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  Tebcm = ott.TmatrixEbcm.simple('sphere', 1.0, ...
      'index_relative', 1.2, 'internal', true);

  Tmie = ott.TmatrixMie(1.0, 'index_relative', 1.2, 'internal', true);

  testCase.verifyThat(Tebcm.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix');

  tol = 1.0e-2;
  Tebcm_data = Tebcm.data ./ max(abs(Tebcm.data(:)));
  Tmie_data = Tmie.data ./ max(abs(Tmie.data(:)));
  testCase.verifyThat(Tebcm_data, IsEqualTo(Tmie_data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance');
    
end