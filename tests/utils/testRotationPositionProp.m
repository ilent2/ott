function tests = testRotationPositionProp
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testNargoutCheck(testCase)

  beam = ott.utils.RotationPositionProp();
  
  testCase.verifyWarning(@() beam.rotateX(5), ...
    'ott:utils:nargoutCheck:no_outputs');
end

function testLocalTranslation(testCase)

  obj = ott.utils.RotationPositionProp();
  obj = obj.rotateY(pi/2);
  obj = obj.translateXyz([0;0;1], 'origin', 'local');
  
  testCase.verifyEqual(obj.position, [1;0;0], 'AbsTol', 1.0e-15);
  testCase.verifyEqual(obj.local2global([0;0;0]), obj.position, 'loc2glob');
  testCase.verifyEqual(obj.global2local(obj.position), [0;0;0], 'glob2loc');

end

function testRotation(testCase)

  obj = ott.utils.RotationPositionProp();
  obj.position = [1;0;0];
  
  trial = obj.rotateY(pi/2, 'origin', 'global');
  testCase.verifyEqual(trial.position, [0;0;-1], 'AbsTol', 1.0e-15, 'global');
  
  trial = obj.rotateY(pi/2, 'origin', 'local');
  testCase.verifyEqual(trial.position, obj.position, 'AbsTol', 1.0e-15, 'local');

end

function testHetrogeneousArray(testCase)

  beam(1) = ott.beam.Empty();
  beam(2) = ott.beam.BscBeam();
  
  tbeam = beam.translateXyz([1;0;0]);
  testCase.verifyEqual([tbeam.position], [1,1;0,0;0,0], 'tx');
  
  tbeam = beam.rotateY(pi/2);
  Ry = ott.utils.roty(pi/2 * 180/pi);
  testCase.verifyEqual([tbeam.rotation], [Ry, Ry], 'Ry');

end
