function tests = testConsistency
  % This file contains tests to verify consistency of Bsc classes
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testEmptyBeams(testCase)
  % Test all beams have an empty method and default constructors have
  % consistent behaviour across all beams.

  % These beams should have an empty method
  % However, their default constructor is non-empty
  partial = {'Gaussian', 'LaguerreGaussian', 'HermiteGaussian', ...
      'InceGaussian', 'PlaneWave', 'Bessel', 'Annular', 'Scattered'};

  for ii = 1:numel(partial)
    beam1 = ott.beam.vswf.(parial{ii});
    beam2 = ott.beam.vswf.(partial{ii}).empty();
    testCase.verifyNotEmpty(beam1, [partial{ii}, ': non-empty']);
    testCase.verifyEmpty(beam2, [partial{ii}, ': empty']);
  end

  % These beams should have an empty method with the same behaviour
  % as the default constructor.
  fully = {'Bsc', 'BesselBasis', 'PlaneBasis', 'LgParaxialBasis'};

  for ii = 1:numel(fully)
    beam1 = ott.beam.vswf.(fully{ii});
    beam2 = ott.beam.vswf.(fully{ii}).empty();
    testCase.verifyEqual(beam1, beam2, [fully{ii}, ': match']);
    testCase.verifyEmpty(beam1, [fully{ii}, ': empty']);
  end
end

