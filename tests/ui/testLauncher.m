function tests = testLauncher
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testLauncherApp(testCase)
  gui = ott.ui.Launcher();
  testCase.addTeardown(@delete,gui);
end