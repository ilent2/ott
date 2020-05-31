function tests = testAnnular
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../../');
end


function testConstruct(testCase)

  angles = [pi/8, pi/4];
  beam = ott.beam.abstract.Annular(angles);

  testCase.verifyEqual(beam.angles, angles);

  h = figure();
  beam.visualiseFarfield();
  close(h);

end

