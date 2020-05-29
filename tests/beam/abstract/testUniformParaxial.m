function tests = testUniformParaxial
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruct(testCase)

  polarisation = [1; 2];
  beam = ott.beam.abstract.UniformParaxial(polarisation);

  testCase.verifyEqual(beam.polarisation, polarisation);

  h = figure();
  beam.visualiseFarfield();
  caxis([0, 10]);
  close(h);

end
