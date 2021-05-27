classdef testDrag < matlab.uitest.TestCase
  
  methods (TestClassSetup)
    function setupPath(testCase)
      addpath('../../');
    end
  end
  
  methods (Test)
    function test_Simple(testCase)

      % Create GUI
      gui = ott.ui.drag.Simple();
      testCase.addTeardown(@delete,gui);

      % Create shape
      ourshape = ott.shape.Sphere();
      assignin('base', 'ourshape', ourshape);
      testCase.addTeardown(@evalin, 'base', 'clear(''ourshape'')');

      % Generate target
      viscosity = 1;
      target = ott.drag.Stokes.FromShape(ourshape, viscosity);

      % Set fields (sphere)
      testCase.addTeardown(@evalin, 'base', 'clear(''output'')');
      testCase.type(gui.VariableName.EditField, 'output');
      testCase.type(gui.ShapeDropDown.DropDown, 'ourshape');
      testCase.type(gui.ViscositySpinner.Spinner, viscosity);

      % Compare output
      output = evalin('base', 'output');
      testCase.verifyEqual(output, target);
    end
  end
end
