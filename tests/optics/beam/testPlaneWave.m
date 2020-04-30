function tests = testPlaneWave
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  direction = [0; 0; 1];
  polarisation = [1; 0; 0];
  beam = ott.optics.beam.PlaneWave(direction, polarisation);
  
  testCase.verifyEqual(beam.direction, direction, 'dir');
  testCase.verifyEqual(beam.polarisation, polarisation, 'pol');
end

function testConstructorArray(testCase)

  direction = [0, 0; 0, 0; 1, -1];
  polarisation = [1; 0; 0];
  beam = ott.optics.beam.PlaneWave(direction, polarisation);
  
  testCase.verifyEqual(beam.direction, direction, 'dir');
  testCase.verifyEqual(beam.polarisation, repmat(polarisation, 1, 2), 'pol');
end

function testVisualise(tsetCase)

  direction = [0, 0; 0, 0; 1, -1];
  polarisation = [1; 0; 0];
  beam = ott.optics.beam.PlaneWave(direction, polarisation);
  
  h = figure();
  beam.visualise();
%   close(h);
  
end

function testFarfield(testCase)

  direction = [0; 0; 1];
  polarisation = [1; 0; 0];
  beam = ott.optics.beam.PlaneWave(direction, polarisation);
  
  h = figure();
  beam.visualiseFarfield();
%   close(h);

  h = figure();
  beam.visualiseFarfield('method', 'delta');
%   close(h);
end
