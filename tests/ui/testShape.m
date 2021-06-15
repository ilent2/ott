classdef testShape < matlab.uitest.TestCase
  
  methods (TestClassSetup)
    function setupPath(testCase)
      addpath('../../');
    end
  end
  
  methods (Test)
    function testCadFile(testCase)
      gui = ott.ui.shape.CadFileLoader();
      testCase.addTeardown(@delete,gui);
      
      fname = '..\shape\cube.stl';
      testCase.type(gui.CadFileSelector.EditField, fname);
    end

    function testSimple(testCase)
      gui = ott.ui.shape.Simple();
      testCase.addTeardown(@delete,gui);
    end

    function testVisualise(testCase)
      gui = ott.ui.shape.Visualise();
      testCase.addTeardown(@delete,gui);
    end
  end
end
