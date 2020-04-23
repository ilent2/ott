function tests = testGaussianParaxial
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  waist = 1.0;
  beam = ott.optics.beam.GaussianParaxial(waist);
  
end

function testVisualise(testCase)

  waist = 1.0;
  beam = ott.optics.beam.GaussianParaxial(waist);
  
  h = figure();
  beam.visualise();
%   close(h);
  
end
