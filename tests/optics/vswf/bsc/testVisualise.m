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

