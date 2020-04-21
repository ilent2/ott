function tests = testEuler
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

force_method = @(position, rotation, time, drag) -[position(:); 0;0;0];

drag = ott.drag.Stokes('forward', 0.1*eye(6));

dynamicsSystem = ott.dynamics.DynamicsSystem.simple(...
  force_method, 'drag', drag, 'brownian_motion', false, ...
  'temperature', 1e-3);

tspan = [0, 1];
position = [1;0;-1];
rotation = eye(3);

simulation = ott.dynamics.MatlabOde(dynamicsSystem, 'solver', @ode23);

[tout, posOut, rotOut, forceOut] = simulation.evaluate(tspan, position, rotation);

end

