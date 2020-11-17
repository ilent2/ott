function tests = testBeam
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testProperties(testCase)

  testCase.verifyEqual(ott.beam.Beam.speed0, 3e8, 'speed0');
  
  % Other properties are tested in testEmpty.m

end

function testVisNearfield(testCase)

  f = figure();
  testCase.addTeardown(@close, f);

  abeam = ott.beam.BscBeam(ott.bsc.Bsc([1; 2; 3], zeros(3, 1)));
  abeam.visNearfield();

  abeam.visNearfield('axis', 'x', 'range', [5, 5]);
  abeam.visNearfield('axis', {[1;0;0], [0;1;0]}, 'range', [-1, 1, -1, 1]);
  abeam.visNearfield('axis', {[1;0;0], [0;1;0], [0;0;1]}, ...
      'range', {linspace(0, 1), linspace(0, 2)});
    
  im = abeam.visNearfield();
  testCase.verifySize(im, [80, 80]);

end

function testVisFarfield(testCase)

  f = figure();
  testCase.addTeardown(@close, f);

  abeam = ott.beam.BscBeam(ott.bsc.Bsc([1; 2; 3], zeros(3, 1)));
  abeam.visFarfield();
    
  im = abeam.visFarfield();
  testCase.verifySize(im, [80, 80]);

end

function testVisFarfieldSlice(testCase)

  f = figure();
  testCase.addTeardown(@close, f);

  abeam = ott.beam.BscBeam(ott.bsc.Bsc([1; 2; 3], zeros(3, 1)));
  abeam.visFarfieldSlice(0.0);

end

function testVisFarfieldSphere(testCase)

  f = figure();
  testCase.addTeardown(@close, f);

  abeam = ott.beam.BscBeam(ott.bsc.Bsc([1; 2; 3], zeros(3, 1)));
  abeam.visFarfieldSphere();

end

function testFarfieldSphereModesMatch(testCase)

  beam1 = ott.beam.BscBeam(ott.bsc.Bsc([1;0;0], zeros(3, 1)));
  beam2 = ott.beam.BscBeam(ott.bsc.Bsc([0;0;1], zeros(3, 1)));
  
  s1 = beam1.visFarfieldSphere();
  s2 = beam2.visFarfieldSphere();
  
  testCase.verifyEqual(s1, s2, 'AbsTol', 1e-16);

end

function testIntensityMoment(testCase)

  beam1 = ott.beam.Gaussian('power', 2);
  [moment, ints] = beam1.intensityMoment();
  emptyBeam = ott.beam.Empty();
  
  testCase.verifyEqual(ints, beam1.power, 'RelTol', 1e-2, 'power');
  testCase.verifyEqual(moment./beam1.speed, ...
    emptyBeam.force(beam1), 'RelTol', 1e-2, 'AbsTol', 1e-23, 'moment');

end

function testForce(testCase)

  beam1 = ott.beam.Gaussian();
  beam2 = ott.beam.BscBeam();
  Fb = beam1.force(beam2);
  
  bsc1 = ott.bsc.Bsc(beam1);
  bsc2 = ott.bsc.Bsc(beam2);
  Fq = bsc1.force(bsc2) * beam1.power ./ beam1.speed;
  
  testCase.verifyEqual(Fb, Fq, 'force');

end

function testTranslateFarfield(testCase)

  beam = ott.beam.Gaussian();
  bsc1 = ott.bsc.Bsc(beam, beam.data.Nmax);
  beam.position = [-0.3;0.2;0.5]*beam.wavelength;
  bsc2 = ott.bsc.Bsc(beam, beam.data.Nmax + 50);
  
  beam1 = ott.beam.BscBeam(bsc1);
  beam2 = ott.beam.BscBeam(bsc2);
  
  % Incoming hemisphere only
  rtp = [ones(1, 5); pi-pi/2*rand(1, 5); 2*pi*rand(1, 5)];
  
  Ertp1 = beam1.efarfield(rtp, 'basis', 'outgoing');
  Ertp2 = beam2.efarfield(rtp, 'basis', 'outgoing');
  
  % Translate Ertp1 to location 2
  Ertp3 = beam1.translateFarfields(Ertp1, -beam.position);
  
  testCase.verifyEqual(Ertp3.vxyz, Ertp2.vxyz, ...
      'AbsTol', 1e-12, 'translated');

end


