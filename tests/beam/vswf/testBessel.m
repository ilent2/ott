function tests = testBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructDefault(testCase)

  beam = ott.beam.vswf.Bessel();

  testCase.verifyEqual(beam.angle, pi/4, 'angle');

end

function testConstructOptional(testCase)

  angle = 0.1;
  beam = ott.beam.vswf.Bessel(angle);

  testCase.verifyEqual(beam.angle, angle, 'angle');

end

