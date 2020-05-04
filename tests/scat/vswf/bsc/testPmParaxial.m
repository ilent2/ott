function tests = bscpmparaxial
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../../');
end

function testConstruct(testCase)

  import ott.optics.vswf.*;
  [xx,yy] = meshgrid(linspace(-1, 1, 128), linspace(-1, 1, 128));
  incident_mode = ott.utils.lgmode(0,0, ...
    sqrt(xx.^2+yy.^2), atan2(yy,xx));

  beam = bsc.PmParaxial(1.0, incident_mode);

end

function testCoefficients(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  import ott.optics.vswf.*;
  tol = 1.0e-8;

  [xx,yy] = meshgrid(linspace(-1, 1, 128), linspace(-1, 1, 128));
  incident_mode = ott.utils.lgmode(0,0, ...
    sqrt(xx.^2+yy.^2), atan2(yy,xx));

  % Create the first beam
  beam1 = bsc.PmParaxial(1.0, incident_mode, ...
    'keep_coefficient_matrix', true, ...
    'invert_coefficient_matrix', false);

  % Create the second beam (reuse data)
  beam2 = bsc.PmParaxial(1.0, incident_mode, 'beamData', beam1);

  % Create the second beam (reuse pass coefficient matrix)
  beam3 = bsc.PmParaxial(1.0, incident_mode, ...
    'beamData', beam1.inv_coefficient_matrix);

  testCase.verifyThat(beam2.getCoefficients(), ...
      IsEqualTo(beam1.getCoefficients(), ...
      'Within', AbsoluteTolerance(tol)), ...
      'beam2 does not match original beam');

  testCase.verifyThat(beam3.getCoefficients(), ...
      IsEqualTo(beam1.getCoefficients(), ...
      'Within', AbsoluteTolerance(tol)), ...
      'beam2 does not match original beam');

end

function testPolarisationNoEffect(testCase)
  % Polarisation should have no effect if input is NxMx2 field
  
  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  import ott.optics.vswf.*;
  tol = 1.0e-8;

  [xx,yy] = meshgrid(linspace(-1, 1, 128), linspace(-1, 1, 128));
  incident_mode = ott.utils.lgmode(0,0, ...
    sqrt(xx.^2+yy.^2), atan2(yy,xx));

  % Create the first beam
  beam1 = bsc.PmParaxial(1.0, incident_mode, ...
    'polarisation', [1, 0]);
  
  % Make incident mode MxNx2
  incident_mode(:, :, 2) = 0;

  % Create the second beam
  beam2 = bsc.PmParaxial(1.0, incident_mode, ...
    'polarisation', [0, 0]);

  testCase.verifyThat(beam2.getCoefficients(), ...
      IsEqualTo(beam1.getCoefficients(), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Beams do not match');
end
  
  
