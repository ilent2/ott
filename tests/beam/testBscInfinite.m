function tests = testBscFinite
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  Nmax = 20;
  dir = [0;0;1];
  pol = [1;0;0];
  bsc = ott.bsc.PlaneWave.FromDirection(Nmax, dir, pol);
  beam = ott.beam.BscInfinite(bsc);

  testCase.assertEqual(beam.data, bsc, 'data not set');

end

function testFieldValues(testCase)

  Nmax = 20;
  dir = [0;0;1];
  pol = [1;0;0];
  bsc = ott.bsc.PlaneWave.FromDirection(Nmax, dir, pol);
  beam = ott.beam.BscInfinite(bsc);

  testCase.verifyWarning(@() beam.efarfield([0;0;0]), ...
    'ott:beam:BscInfinite:farfield_is_finite', 'efarfield');
  testCase.verifyWarning(@() beam.hfarfield([0;0;0]), ...
    'ott:beam:BscInfinite:farfield_is_finite', 'hfarfield');
  testCase.verifyWarning(@() beam.ehfarfield([0;0;0]), ...
    'ott:beam:BscInfinite:farfield_is_finite', 'ehfarfield');
  
  large = [1;1;1];
  testCase.verifyError(@() beam.efield(large), ...
    'ott:beam:BscBeam:recalculate_not_implemented', 'efield');
  testCase.verifyError(@() beam.hfield(large), ...
    'ott:beam:BscBeam:recalculate_not_implemented', 'hfield');
  testCase.verifyError(@() beam.ehfield(large), ...
    'ott:beam:BscBeam:recalculate_not_implemented', 'ehfield');
  
  testCase.verifyError(@() beam.efieldRtp(large), ...
    'ott:beam:BscBeam:recalculate_not_implemented', 'efieldRtp');
  testCase.verifyError(@() beam.hfieldRtp(large), ...
    'ott:beam:BscBeam:recalculate_not_implemented', 'hfieldRtp');
  testCase.verifyError(@() beam.ehfieldRtp(large), ...
    'ott:beam:BscBeam:recalculate_not_implemented', 'ehfieldRtp');

end


