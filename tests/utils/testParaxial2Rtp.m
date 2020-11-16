function tests = testParaxial2Rtp
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testAngleRanges(testCase)

  [x, y] = meshgrid(linspace(-1, 1), linspace(-1, 1));
  xy = [x(:), y(:)].';
  mapping = 'sin';
  direction = 'pos';
  rtp = ott.utils.paraxial2rtp(xy, mapping, direction);
  
  filt = ~isnan(rtp(2, :));
  
  testCase.verifyGreaterThanOrEqual(rtp(1, filt), 0, 'r');
  testCase.verifyGreaterThanOrEqual(rtp(2, filt), 0, 't');
  testCase.verifyLessThanOrEqual(rtp(2, filt), pi, 't pi');
  testCase.verifyGreaterThanOrEqual(rtp(3, filt), 0, 'p');
  testCase.verifyLessThan(rtp(3, filt), 2*pi, 'p 2pi');

end
