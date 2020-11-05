function tests = testWBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testCast(testCase)

  beam = ott.beam.Bessel();
  testCase.assertInstanceOf(beam, 'ott.beam.BscWBessel');
  
  Nmax = 5;
  bsc = ott.bsc.Bsc(beam, Nmax);

end

