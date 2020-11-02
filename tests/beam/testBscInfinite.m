function tests = testBscFinite
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  bsc = ott.bsc.Bsc(ones(8, 1), 0);
  beam = ott.beam.BscInfinite(bsc);

  testCase.assertEqual(beam.data, bsc, 'data not set');

  % Check large get data
  testCase.verifyError(beam.getData(3), ...
      'ott:beam:BscInfinite:Nmax_outside_range');

  % Check large translation
  testCase.verifyErrorFree(beam.getData(2));
  beam.position = [1;0;0];
  testCase.verifyError(beam.getData(2), ...
      'ott:beam:BscInfinite:Nmax_outside_range');

end

