function tests = testVisualise
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testVisualiseFunction(testCase)

  beam = ott.beam.vswf.PmGauss();
  beam.basis = 'regular';
  
  h = figure();
  beam.visualise();
  close(h);

end

function testVisualiseFarfield(testCase)

  beam = ott.beam.vswf.PmGauss();
  beam.basis = 'incoming';
  
  h = figure();
  beam.visualiseFarfield('dir', 'neg');
  close(h);

end

function testVisualiseFarfieldSlice(testCase)

  beam = ott.beam.vswf.PmGauss();
  beam.basis = 'incoming';
  phi = 0.0;
  
  h = figure();
  beam.visualiseFarfieldSlice(phi);
  close(h);

end

function testVisualiseFarfieldSphere(testCase)

  beam = ott.beam.vswf.PmGauss();
  beam.basis = 'incoming';
  
  h = figure();
  beam.visualiseFarfieldSphere();
  close(h);

end

function testVisualiseArray(testCase)

  beam = ott.beam.vswf.PmGauss();
  beam = [beam, beam];
  
  im = beam.visualise();
  testCase.verifyEqual(size(im), [1, 2], 'default');
  testCase.verifyEqual(size(im{1}), [80, 80], 'default cell');

  beam.array_type = 'coherent';
  im = beam.visualise();
  testCase.verifyEqual(size(im), [80, 80], 'coherent');
  
  beam.array_type = 'incoherent';
  im = beam.visualise();
  testCase.verifyEqual(size(im), [80, 80], 'incoherent');
  
  beam = [beam, beam(1)];
  beam = ott.beam.Array('array', [1, 2], beam, beam);
  
  im = beam.visualise();
  testCase.verifyEqual(size(im), [1, 2], 'matrix');
  testCase.verifyEqual(size(im{1}), [80, 80], 'matrix cell');
end

