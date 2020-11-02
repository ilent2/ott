function tests = testBscFinite
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  bsc = ott.bsc.Bsc([1;0;0], 0);
  beam = ott.beam.BscFinite(bsc);

  beam.position = [1;0;0]*1e-6;
  Nmax = 10;
  tbeam = beam.getData(Nmax);

  tbsc = bsc.translateXyz(beam.position, 'Nmax', Nmax);

  testCase.verifyEqual(tbeam, tbsc, 'bsc data doesnt match');

end

