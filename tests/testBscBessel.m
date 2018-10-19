function tests = bscbessel
  tests = functiontests(localfunctions);
end

function testConstruct(testCase)

  addpath('../');
  beam = ott.BscBessel(12, [0.0, pi/4]);

  import matlab.unittest.constraints.IsEqualTo;
  testCase.verifyThat(beam.Nbeams, IsEqualTo(2), ...
    'Incorrect number of beams stored');

end

