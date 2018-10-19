function tests = bscpmparaxial
  tests = functiontests(localfunctions);
end

function testConstruct(testCase)

  addpath('../');
  
  [xx,yy] = meshgrid(linspace(-1, 1, 512), linspace(-1, 1, 512));
  incident_mode = ott.utils.lgmode(0,0, ...
    sqrt(xx.^2+yy.^2), atan2(yy,xx));
  
  beam = ott.BscPmParaxial(1.0, incident_mode);

end

