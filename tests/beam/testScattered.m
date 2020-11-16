function tests = testScattered
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstructDefault(testCase)

  beam = ott.beam.Scattered();
  testCase.verifyEqual(beam.scattered, [], 'scattered');
  testCase.verifyEqual(beam.incident, ott.beam.Empty(), 'incident');
  testCase.verifyEqual(beam.particle, [], 'particle');
  testCase.verifyEqual(beam.internal, [], 'internal');
  
  beam.defaultVisRangeInternal();

end

function testFieldFunctionsEmpty(testCase)

  beam = ott.beam.Scattered();
  rtp = [zeros(3, 1), [10;0;0]];
  xyz = ott.utils.rtp2xyz(rtp);
  
  target = zeros(3, 2);
  testCase.verifyEqual(beam.efieldRtp(rtp).vxyz, target, 'er no shape');
  testCase.verifyEqual(beam.efield(xyz).vxyz, target, 'e no shape');
  testCase.verifyEqual(beam.hfieldRtp(rtp).vxyz, target, 'hr no shape');
  testCase.verifyEqual(beam.hfield(xyz).vxyz, target, 'h no shape');
  testCase.verifyEqual(beam.ehfieldRtp(rtp).vxyz, target, 'ehr no shape');
  testCase.verifyEqual(beam.ehfield(xyz).vxyz, target, 'eh no shape');
  
  beam.shape = ott.shape.Sphere(1);
  target = [nan(3, 1), zeros(3, 1)];
  testCase.verifyEqual(beam.efieldRtp(rtp).vxyz, target, 'er shape');
  testCase.verifyEqual(beam.efield(xyz).vxyz, target, 'e shape');
  testCase.verifyEqual(beam.hfieldRtp(rtp).vxyz, target, 'hr shape');
  testCase.verifyEqual(beam.hfield(xyz).vxyz, target, 'h shape');
  testCase.verifyEqual(beam.ehfieldRtp(rtp).vxyz, target, 'ehr shape');
  testCase.verifyEqual(beam.ehfield(xyz).vxyz, target, 'eh shape');
  
end

function testFarfieldEmpty(testCase)

  beam = ott.beam.Scattered();
  rtp = [10;0;0];
  xy = [0;0];

  % Test coverage
  beam.efarfield(rtp);
  beam.hfarfield(rtp);
  beam.ehfarfield(rtp);
  beam.hparaxial(xy);
  beam.eparaxial(xy);
  beam.ehparaxial(xy);

end

function testFieldPartialEmpty(testCase)

  beam1 = ott.beam.PlaneWave('Nmax', 10);
  beam2 = beam1.rotateY(pi/2);
  beam3 = beam1.rotateX(pi/2);
  
  xyz = [1;2;3]*beam1.wavelength;
  E1 = beam1.efield(xyz).vxyz;
  E2 = beam2.efield(xyz).vxyz;
  E3 = beam3.efield(xyz).vxyz;
  
  scat = ott.beam.Scattered('incident', beam1, 'scattered', beam2);
  E = scat.efield(xyz, 'type', 'total').vxyz;
  testCase.verifyEqual(E, 2*E2 + E1, 'total = I + 2*S');
  
  scat = ott.beam.Scattered('incident', beam1, 'internal', beam3);
  E = scat.efield(xyz, 'type', 'total').vxyz;
  testCase.verifyEqual(E, E3, 'internal only');
  
end

function testForceTorque(testCase)

  beam = ott.beam.Scattered('scattered', ott.beam.Empty());
  
  testCase.verifyEqual(beam.force(), zeros(3, 1), 'force');
  beam.torque();
  beam.spin();
  beam.forcetorque();
  
end

function testContinuousFields(testCase)

  radius = 1.5;
  index_rel = 1.5;
  index_med = 1.5;
  
  [texternal, tinternal] = ott.tmatrix.Mie(radius, 'index_relative', index_rel);
  
  pexternal = ott.particle.Fixed('tmatrix', texternal);
  pinternal = ott.particle.Fixed('tinternal', tinternal);
  
  ibeam = ott.beam.PlaneWave('Nmax', texternal.Nmax(2)+5, 'index_medium', index_med);
  sbeam = ibeam.scatter(pexternal);
  nbeam = ibeam.scatter(pinternal);
  
  % Locations to calculate fields (boudnary)
  N = [[1;0;0], [0;0;1]];
  xyz = N*radius*ibeam.wavelength;
  
  % Calculate fields
  [Es, Hs] = sbeam.ehfield(xyz);
  [En, Hn] = nbeam.ehfield(xyz);
  
  n2 = index_rel^2;
  m2 = index_med^2;
  
  testCase.verifyEqual(cross(N, En.vxyz), cross(N, Es.vxyz), ...
    'RelTol', 1e-3, ...
    'tangential electric field (E) should be continuous');
  
  testCase.verifyEqual(dot(N, n2*En.vxyz), dot(N, Es.vxyz), ...
    'RelTol', 1.5e-4, ...
    'perpendicular displacement field (D) should be continuous');
  
  testCase.verifyEqual(cross(N, Hn.vxyz), cross(N, Hs.vxyz), ...
    'RelTol', 1e-3, ...
    'tangential magnetic field (B) should be continuous');
  
  % No refractive index term (mu = mu0 in both)
  testCase.verifyEqual(dot(N, Hn.vxyz), dot(N, Hs.vxyz), ...
    'RelTol', 2e-4, ...
    'perpendicular magnetic field (H) should be continuous');
end

function testScatteringPreservesPower(testCase)

  beam = ott.beam.Gaussian();
  [~, I0] = beam.intensityMoment();
  
  particle = ott.particle.Fixed.FromShape(...
      ott.shape.Sphere(beam.wavelength*1.1), 'index_relative', 1.1, ...
      'wavelength0', beam.wavelength0);
  sbeam = beam.scatter(particle);
  [~, I1] = sbeam.intensityMoment();
  
  testCase.verifyEqual(I1, I0, 'RelTol', 1e-3, ...
    'beam power should be preserved with non-absorbing particle');

end

function testFarfieldTranslationsIncident(testCase)

  beam1 = ott.beam.Gaussian();
  beam2 = ott.beam.Empty();
  beam = ott.beam.Scattered('incident', beam1, 'scattered', beam2);
  
  rtp = [ones(1, 5); mod(abs(randn(2, 5)), pi)];
  
  Ertp1 = beam1.efarfield(rtp);
  Ertp2 = beam.efarfield(rtp);
  testCase.assertEqual(Ertp1.vxyz, Ertp2.vxyz, 'original position');
  
  % Tranlsate incident and beam
  xyz = randn(3, 1)*beam.wavelength;
  beam.position = xyz;
  beam1.position = xyz;
  
  Ertp1 = beam1.efarfield(rtp);
  Ertp2 = beam.efarfield(rtp);
  testCase.verifyEqual(Ertp1.vxyz, Ertp2.vxyz, 'scat.position');
  
  xyz2 = randn(3, 1)*beam.wavelength;
  beam.incident.position = xyz2;
  beam1.position = xyz + xyz2;
  
  Ertp1 = beam1.efarfield(rtp);
  Ertp2 = beam.efarfield(rtp);
  testCase.verifyEqual(Ertp1.vxyz, Ertp2.vxyz, ...
    'AbsTol', 2e-14, 'scat.incident');
  
end

function testFarfieldTranslationsScattered(testCase)

  beam1 = ott.beam.Gaussian();
  beam2 = ott.beam.Empty();
  beam = ott.beam.Scattered('incident', beam2, 'scattered', beam1);
  beam1 = beam1 * 2;
  
  rtp = [ones(1, 5); mod(abs(randn(2, 5)), pi)];
  
  Ertp1 = beam1.efarfield(rtp);
  Ertp2 = beam.efarfield(rtp);
  testCase.assertEqual(Ertp1.vxyz, Ertp2.vxyz, 'original position');
  
  % Tranlsate incident and beam
  xyz = randn(3, 1)*beam.wavelength;
  beam.position = xyz;
  beam1.position = xyz;
  
  Ertp1 = beam1.efarfield(rtp);
  Ertp2 = beam.efarfield(rtp);
  testCase.verifyEqual(Ertp1.vxyz, Ertp2.vxyz, 'scat.position');
  
  xyz2 = randn(3, 1)*beam.wavelength;
  beam.scattered.position = xyz2;
  beam1.position = xyz + xyz2;
  
  Ertp1 = beam1.efarfield(rtp);
  Ertp2 = beam.efarfield(rtp);
  testCase.verifyEqual(Ertp1.vxyz, Ertp2.vxyz, ...
    'AbsTol', 2e-14, 'scat.scattered');
  
end
