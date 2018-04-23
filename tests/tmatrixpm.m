function tests = tmatrixpm
  tests = functiontests(localfunctions);
end

function testConstruct(testCase)

  addpath('../');
  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  Tpm = ott.TmatrixPm.simple('sphere', 1.0, 'index_relative', 1.2);

  Tmie = ott.TmatrixMie(1.0, 'index_relative', 1.2);

  testCase.verifyThat(Tpm.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix');

  tol = 1.0e-2;
  Tpm_data = Tpm.data ./ max(abs(Tpm.data(:)));
  Tmie_data = Tmie.data ./ max(abs(Tmie.data(:)));
  testCase.verifyThat(Tpm_data, IsEqualTo(Tmie_data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance');

end

