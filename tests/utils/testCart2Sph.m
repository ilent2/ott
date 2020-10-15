function tests = testCart2Sph
  % Test Cartesian and Spherical coordinate transform functions
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testCoordinateTransform(testCase)

  oxyz = [1; 2; 3];
  rtp = ott.utils.xyz2rtp(oxyz);
  xyz = ott.utils.rtp2xyz(rtp);

  testCase.verifyEqual(xyz, oxyz, ...
    'AbsTol', 1e-15, 'Forward-back transform not equal');
end

function testShapePreserved(testCase)

  xyz = rand(3, 4);

  [r, ~, ~] = ott.utils.xyz2rtp(xyz(1, :), xyz(2, :), xyz(3, :));
  testCase.verifyEqual(size(r), [1, 4]);

  [r, ~, ~] = ott.utils.xyz2rtp(xyz(1, :).', xyz(2, :).', xyz(3, :).');
  testCase.verifyEqual(size(r), [4, 1]);

  [r, ~, ~] = ott.utils.rtp2xyz(xyz(1, :), xyz(2, :), xyz(3, :));
  testCase.verifyEqual(size(r), [1, 4]);

  [r, ~, ~] = ott.utils.rtp2xyz(xyz(1, :).', xyz(2, :).', xyz(3, :).');
  testCase.verifyEqual(size(r), [4, 1]);

end

function testVectorTransform(testCase)

  oxyz = [1; 2; 3]; ovxyz = [0; 0; 1];
  [rtp_vec, rtp_pos] = ott.utils.xyzv2rtpv(ovxyz, oxyz);
  [vxyz, xyz] = ott.utils.rtpv2xyzv(rtp_vec, rtp_pos);

  testCase.verifyEqual(vxyz, ovxyz, ...
    'AbsTol', 1e-15, 'Forward-back transform not equal (vector)');
  testCase.verifyEqual(xyz, oxyz, ...
    'AbsTol', 1e-15, 'Forward-back transform not equal (position)');

end

function testVectorTransformManyInputs(testCase)

  oxyz = [1; 2; 3]; ovxyz = [0; 0; 1];
  [rtp_vec, rtp_pos] = ott.utils.xyzv2rtpv(ovxyz(1), ovxyz(2), ...
    ovxyz(3), oxyz(1), oxyz(2), oxyz(3));
  [vxyz, xyz] = ott.utils.rtpv2xyzv(rtp_vec(1), rtp_vec(2), ...
    rtp_vec(3), rtp_pos(1), rtp_pos(2), rtp_pos(3));

  testCase.verifyEqual(vxyz, ovxyz, ...
    'AbsTol', 1e-15, 'Forward-back transform not equal (vector)');
  testCase.verifyEqual(xyz, oxyz, ...
    'AbsTol', 1e-15, 'Forward-back transform not equal (position)');

end
