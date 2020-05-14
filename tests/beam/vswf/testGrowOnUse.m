function tests = testGrowOnUse
  % Tests for GrowOnUse and classes inheriting from this class
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testGrowOnUse(testCase)

  % List of grow-on-use beams
  beam_classes = {'Gaussian', 'LaguerreGaussian', 'HermiteGaussian', ...
      'InceGaussian', 'PlaneWave', 'Bessel', 'Annular'};

  for ii = 1:numel(beam_classes)

    beam = ott.beam.vswf.(beam_classes{ii})();
    testCase.verifyEqual(beam.Nmax, 0, [beam_classes{ii}, ': Nmax']);
    testCase.verifyNotEmpty(beam, [beam_classes{ii}, ': empty']);

    beam2 = ott.beam.vswf.(beam_classes{ii})('Nmax', 10);
    beam.Nmax = beam2.Nmax;
    testCase.verifyEqual(beam, beam2, ['new beam: ', beam_classes{ii}]);
  end
end

