classdef testTools < matlab.uitest.TestCase
  
  methods (TestClassSetup)
    function setupPath(testCase)
      addpath('../../');
    end
  end
  
  methods (Test)
    function testForcePosition(testCase)
      gui = ott.ui.tools.ForcePosition();
      testCase.addTeardown(@delete,gui);
      
      beam = ott.beam.Gaussian();
      assignin('base', 'beam', beam);
      testCase.addTeardown(@evalin, 'base', 'clear(''beam'')');
      
      tmatrix = ott.tmatrix.Mie(1, 1.2);
      assignin('base', 'tmatrix', tmatrix);
      testCase.addTeardown(@evalin, 'base', 'clear(''tmatrix'')');
      
      testCase.type(gui.BeamDropDown.DropDown, 'beam');
      testCase.type(gui.TmatrixDropDown.DropDown, 'tmatrix');
      testCase.type(gui.ResolutionSpinner.Spinner, 10);
      testCase.press(gui.UpdateButton.Button);
    end

    function testPowerSpectrum(testCase)
      gui = ott.ui.tools.PowerSpectrum();
      testCase.addTeardown(@delete,gui);
      
      testCase.type(gui.XInputDropDown.DropDown, 'linspace(0,100,1000)');
      testCase.type(gui.YInputDropDown.DropDown, 'randn(1,1000) + sin(10*linspace(0,100,1000))');
    end
  end
end