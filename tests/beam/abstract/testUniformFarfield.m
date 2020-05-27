function tests = testUniformFarfield
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruct(testCase)

  polarisation = [1; 2];
  beam = ott.beam.abstract.UniformFarfield(polarisation);

  testCase.verifyEqual(beam.polarisation, polarisation);

  h = figure();
  beam.visualiseFarfield();
  close(h);

end
