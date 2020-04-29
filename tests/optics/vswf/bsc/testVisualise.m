function tests = testVisualise
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../../');
end

function testVisualiseFunction(testCase)

  beam = ott.optics.vswf.bsc.PmGauss();
  beam.basis = 'regular';
  
  h = figure();
  beam.visualise();
  close(h);

end

function testVisualiseFarfield(testCase)

  beam = ott.optics.vswf.bsc.PmGauss();
  beam.basis = 'incoming';
  
  h = figure();
  beam.visualiseFarfield('dir', 'neg');
  close(h);

end

function testVisualiseFarfieldSlice(testCase)

  beam = ott.optics.vswf.bsc.PmGauss();
  beam.basis = 'incoming';
  phi = 0.0;
  
  h = figure();
  beam.visualiseFarfieldSlice(phi);
  close(h);

end

function testVisualiseFarfieldSphere(testCase)

  beam = ott.optics.vswf.bsc.PmGauss();
  beam.basis = 'incoming';
  
  h = figure();
  beam.visualiseFarfieldSphere();
  close(h);

end

function testVisualiseArray(testCase)

  beam = ott.optics.vswf.bsc.PmGauss();
  beam(2) = beam;
  
  im = beam.visualise();
  testCase.verifyEqual(size(im), [80, 80, 2], 'default');

  im = beam.visualise('combine', 'coherent');
  testCase.verifyEqual(size(im), [80, 80], 'coherent');
  
  im = beam.visualise('combine', 'incoherent');
  testCase.verifyEqual(size(im), [80, 80], 'incoherent');
  
  beam(3) = beam(1);
  beam = [beam; beam];
  
  im = beam.visualise();
  testCase.verifyEqual(size(im), [80, 80, 2, 3], 'matrix');
end

function testVisualiseAppend(testCase)

  beam = ott.optics.vswf.bsc.PmGauss();
  beam = beam.append(beam);
  
  im = beam.visualise();
  testCase.verifyEqual(size(im), [80, 80, 2], 'default');

  im = beam.visualise('combine', 'coherent');
  testCase.verifyEqual(size(im), [80, 80], 'coherent');
  
  im = beam.visualise('combine', 'incoherent');
  testCase.verifyEqual(size(im), [80, 80], 'incoherent');
  
  beam = [beam, beam, beam];
  beam = [beam; beam];
  
  im = beam.visualise();
  testCase.verifyEqual(size(im), [80, 80, 2, 2, 3], 'matrix');
end
