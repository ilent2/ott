function tests = bscpmgauss
  tests = functiontests(localfunctions);
end

function testConstruct(testCase)

  addpath('../');
  beam = ott.BscPmGauss();

  import matlab.unittest.constraints.IsEqualTo;
  testCase.verifyThat(beam.beams, IsEqualTo(1), ...
    'Incorrect number of beams stored');

  testCase.verifyThat(beam.mode, IsEqualTo([0 0]), ...
    'Incorrect number of beams stored');

end

