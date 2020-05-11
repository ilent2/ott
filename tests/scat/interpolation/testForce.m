function tests = testForce
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruction(testCase)

  position = randn(3, 100);
  force = exp(-position(3, :).^2) .* [0;0;1];

  part = ott.scat.interp.Force(position, force);

  pred_position = randn(1, 20) .* [0;0;1];
  pred_force = part.force('position', pred_position);
  
  testCase.verifyEqual(size(pred_force), size(pred_position));
  
%   figure();
%   plot(position(3, :), force(3, :), '.', pred_position(3, :), pred_force(3, :), '*');

end

