classdef testParticle < matlab.uitest.TestCase
  
  methods (TestClassSetup)
    function setupPath(testCase)
      addpath('../../');
    end
  end
  
  methods (Test)
    function testFixed(testCase)
      gui = ott.ui.particle.Fixed();
      testCase.addTeardown(@delete,gui);
    end

    function testFromShape(testCase)
      gui = ott.ui.particle.FromShape();
      testCase.addTeardown(@delete,gui);

      % Parameters
      viscosity = 1;
      mass = 2;
      radius = 3;
      wavelength = 1.5;
      irPart = 1.5;
      irMed = 1.2;
      
      % Create shape
      ourshape = ott.shape.Sphere(radius);
      assignin('base', 'ourshape', ourshape);
      testCase.addTeardown(@evalin, 'base', 'clear(''ourshape'')');
      
      % Generate target
      target = ott.particle.Fixed.FromShape(ourshape, ...
        'viscosity', viscosity, 'mass', mass, 'wavelength0', wavelength, ...
        'index_medium', irMed, 'index_particle', irPart);
      
      % Set fields
      testCase.addTeardown(@evalin, 'base', 'clear(''output'')');
      testCase.type(gui.VariableName.EditField, 'output');
      testCase.type(gui.ShapeDropDown.DropDown, 'ourshape');
      testCase.type(gui.MediumIrSpinner.Spinner, irMed);
      testCase.type(gui.ParticleIrSpinner.Spinner, irPart);
      testCase.type(gui.ViscositySpinner.Spinner, viscosity);
      testCase.type(gui.MassSpinner.Spinner, mass);
      testCase.type(gui.WavelengthSpinner.Spinner, wavelength);
     
      % Compare output
      testCase.press(gui.UpdateButton.Button);
      output = evalin('base', 'output');
      testCase.verifyEqual(output, target);
      
    end
  end
end

