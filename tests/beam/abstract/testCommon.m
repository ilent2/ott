function tests = testCommon
  % Tests to run on multiple abstract beams

  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testCastBeam(testCase)
  % Check that the Beam casts are implemented
  % Doesn't check the values of the output, this should be done in
  % the individual beam test files.

  beams = {'Annular', 'Bessel', 'Dipole', 'Empty', 'FarfieldMasked', ...
      'FocussedTopHat', 'Gaussian', 'HermiteGaussian', 'InceGaussian', ...
      'LaguerreGaussian', 'NearfieldMasked', 'ParaxialMasked', ...
      'PlaneWave', 'Ray', 'TopHat', 'Array', 'Incoherent', 'Coherent', ...
      'Scattered'};

  for ii = 1:length(beams)
    meta = meta.class.fromName(beams{ii});
    methods = {meta.MethodList.Name};
    testCase.verifyThat(any(strcmpi('ott.beam.Beam', methods)), ...
        beams{ii});
  end
end

function testCastVswfBsc(testCase)
  % Check that the Bsc casts are implemented
  % Doesn't check the values of the output, this should be done in
  % the individual beam test files.

  beams = {'Annular', 'Bessel', 'Dipole', 'Empty', 'FarfieldMasked', ...
      'FocussedTopHat', 'Gaussian', 'HermiteGaussian', 'InceGaussian', ...
      'LaguerreGaussian', 'NearfieldMasked', 'ParaxialMasked', ...
      'PlaneWave', 'Ray', 'TopHat', 'Array', 'Incoherent', 'Coherent', ...
      'Scattered'};

  for ii = 1:length(beams)
    meta = meta.class.fromName(beams{ii});
    methods = {meta.MethodList.Name};
    testCase.verifyThat(any(strcmpi('ott.beam.vswf.Bsc', methods)), ...
        beams{ii});
  end
end

function testCastParaxial(testCase)
  % Check that the Paraxial casts are implemented
  % Doesn't check the values of the output, this should be done in
  % the individual beam test files.

  beams = {'Gaussian', 'HermiteGaussian', 'InceGaussian', 'LaguerreGaussian'};

  for ii = 1:length(beams)
    meta = meta.class.fromName(beams{ii});
    methods = {meta.MethodList.Name};
    testCase.verifyThat(any(strcmpi('ott.beam.paraxial.Paraxial', methods)), ...
        beams{ii});
  end
end
