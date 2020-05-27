function tests = testBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
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

