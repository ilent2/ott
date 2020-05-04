function tests = testGaussianDavis5
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  waist = 1.0;
  beam = ott.optics.beam.GaussianDavis5(waist);
  
end

function testVisualise(testCase)

  waist = 0.5;
  beam = ott.optics.beam.GaussianDavis5(waist);
  
  h = figure();
  beam.visualise();
  axis image;
  close(h);
  
end
