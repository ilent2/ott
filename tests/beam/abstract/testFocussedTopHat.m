function tests = testBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruct(testCase)

  angle = pi/4;
  beam = ott.beam.abstract.FocussedTopHat(angle);

  testCase.verifyEqual(beam.angle, angle);

  h = figure();
  beam.visualiseFarfield();
  close(h);

end
