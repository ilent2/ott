function tests = testPlaneWave
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  polarisation = [1, 0];
  Nmax = 2;
  beam = ott.beam.PlaneWave(polarisation, 'Nmax', Nmax);

  testCase.verifyEqual(beam.polarisation, polarisation, 'pol');
  
  testCase.assertClass(beam.data, 'ott.bsc.PlaneWave', 'class data');
  testCase.verifyEqual(beam.Nmax, Nmax, 'Nmax');
  testCase.verifyEqual(beam.data(1).direction, [0;0;1], 'direction');
  
  % Check Nmax grows
  newNmax = 13;
  bsc = ott.bsc.Bsc(beam,newNmax);
  testCase.verifyEqual(bsc.Nmax, newNmax);
  
  % Check set Nmax
  beam.Nmax = 14;
  testCase.verifyEqual(beam.Nmax, 14, 'set nmax');

end

function testConstructZeroNmax(testCase)

  polarisation = [1, 0];
  Nmax = 0;
  beam = ott.beam.PlaneWave(polarisation, 'Nmax', Nmax);
  
  bsc = ott.bsc.Bsc(beam);
  testCase.verifyEqual(bsc.nbeams, 1, 'num beams');
end

