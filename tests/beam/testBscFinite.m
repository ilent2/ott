function tests = testBscFinite
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  index = 1.33;
  omega = 2*pi*3e8 ./ 532e-9;
  bsc = ott.bsc.Bsc([1;2;3], [4;5;6]);
  beam = ott.beam.BscFinite(bsc, 'index_medium', index, 'omega', omega);
  
  testCase.verifyEqual(beam.data, bsc, 'bsc');
  testCase.verifyEqual(beam.index_medium, index, 'index');
  testCase.verifyEqual(beam.omega, omega, 'omega');
end

function testGetData(testCase)

  bsc = ott.bsc.Bsc([1;2;3;4;5;6;7;8], zeros(8, 1));
  beam = ott.beam.BscFinite(bsc);
  
  % Original beam
  trial = ott.bsc.Bsc(beam);
  testCase.verifyEqual(trial, bsc, 'original');
  
  % Grow Nmax (should be ignored)
  trial = ott.bsc.Bsc(beam, 10);
  testCase.verifyEqual(trial, bsc, 'grow');
  
  % Translate and grow
  Nmax = 10;
  tbsc = bsc.translateXyz([0.001;0;0], 'Nmax', Nmax);
  tbeam = beam.translateXyz([0.001;0;0]*beam.wavelength);
  trial = ott.bsc.Bsc(tbeam, Nmax);
  testCase.verifyEqual(trial, tbsc, 'translate');
  
end

