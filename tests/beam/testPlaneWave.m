function tests = testPlaneWave
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstructor(testCase)

  direction = [0; 0; 1];
  polarisation = [1; 0; 0];
  beam = ott.beam.PlaneWave('direction', direction, ...
    'polarisation', polarisation);
  
  testCase.verifyEqual(beam.direction, direction, 'dir');
  testCase.verifyEqual(beam.polarisation, polarisation, 'pol');
end

function testConstructorArray(testCase)

  direction = [0, 0; 0, 0; 1, -1];
  polarisation = [1; 0; 0];
  beam = ott.beam.PlaneWave('direction', direction, ...
    'polarisation', polarisation);
  
  testCase.verifyEqual(beam.direction, direction, 'dir');
  testCase.verifyEqual(beam.polarisation, repmat(polarisation, 1, 2), 'pol');
end

function testFromFarfield(testCase)
  beam = ott.beam.paraxial.Gaussian(1.0);
  beam = ott.beam.PlaneWave.FromFarfield(beam);
end

function testFromNearfield(testCase)
  beam = ott.beam.paraxial.Gaussian(1.0);
  beam = ott.beam.PlaneWave.FromNearfield(beam);
end

function testFromParaxial(testCase)
  beam = ott.beam.paraxial.Gaussian(1.0);
  beam = ott.beam.PlaneWave.FromParaxial(beam);
end

function testVisualise(tsetCase)

  direction = [0, 0; 0, 0; 1, -1];
  polarisation = [1; 0; 0];
  beam = ott.beam.PlaneWave('direction', direction, ...
    'polarisation', polarisation);
  
  h = figure();
  beam.visualise();
  close(h);
  
end

function testFarfield(testCase)

%   gbeam = ott.optics.beam.GaussianParaxial(1.0);
%   beam = ott.optics.beam.PlaneWave.FromFarfield(gbeam);

  direction = [0; 0; 1];
  polarisation = [1; 0; 0];
  beam = ott.beam.PlaneWave('direction', direction, ...
    'polarisation', polarisation);
  
  h = figure();
  beam.visualiseFarfieldSphere();
  close(h);

  h = figure();
  beam.visualiseFarfieldSphere('method', 'delta');
  close(h);
end
