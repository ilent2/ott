% Generate the figures for the Spherical Wave Representation section

%% different kinds of beams

figure();

subplot(1, 3, 1);
a = [0, 0, 0, 0, 1, 0, 0, 0];  % [3 dipole, 5 quadrapole]
b = 1i * a;
basis = 'incoming';
type = 'incident';
beam = ott.Bsc(a, b, basis, type);
beam.visualiseFarfieldSphere('type', '3dpolar');
title('a');

subplot(1, 3, 2);
beam = ott.BscPlane(0, 0, 'Nmax', 12, 'wavelength0', 1064e-9, ...
  'index_medium', 1.0);
x = linspace(-3, 3, 80)*1e-6;
y = x;
beam.visualise('axis', 'x', 'range', {x, y});
drawNmax([0,0], beam.Nmax, 'r', ':');
title('b');

subplot(1, 3, 3);
beam = ott.BscPmGauss('NA', 0.98, 'polarisation', [1, 1i], ...
  'wavelength0', 1064e-9, 'index_medium', 1.0);
beam.Nmax = 10;
x = linspace(-2, 2, 80)*1e-6;
y = x;
beam.visualise('axis', 'x', 'range', {x, y});
drawNmax([0,0], beam.Nmax, 'r', ':');
hold on;
ka = ott.utils.nmax2ka(beam.Nmax);
a = ka ./ (2*pi/1064e-9);
h = plot([0, 0], [-a, a], 'w');
h.LineWidth = 2;
hold off;
title('c');

%% Plane wave translations

beam = ott.BscPlane(0, 0, 'Nmax', 12, 'wavelength0', 1064e-9, ...
  'index_medium', 1.0);

x = linspace(-3, 3, 80)*1e-6;
y = x;

figure();
subplot(2, 3, 1);
beam.visualise('axis', 'x', 'range', {x, y});
drawNmax([0,0], beam.Nmax, 'r', ':');
drawNmax([0,0.5e-6], beam.Nmax/2, 'b', ':');
drawNmax([0,1.9e-6], beam.Nmax/2, 'w', ':');
% caxr = caxis();
caxr = [0, 2];
caxis(caxr);
title('a');

subplot(2, 3, 2);
tbeam = beam.translateXyz([-0.5e-6; 0; 0], 'Nmax', beam.Nmax/2);
tbeam.visualise('axis', 'x', 'range', {x, y});
drawNmax([0,0], tbeam.Nmax, 'b', ':');
drawNmax([0,-0.5e-6], beam.Nmax, 'r', ':');
caxis(caxr);
title('b');

subplot(2, 3, 3);
tbeam = beam.translateXyz([-1.9e-6; 0; 0], 'Nmax', beam.Nmax/2);
tbeam.visualise('axis', 'x', 'range', {x, y});
drawNmax([0,0], tbeam.Nmax, 'w', ':');
drawNmax([0,-1.9e-6], beam.Nmax, 'r', ':');
caxis(caxr);
title('c');

%% Gaussian translations

beam = ott.BscPmGauss('NA', 0.98, 'polarisation', [1, 1i], ...
  'wavelength0', 1064e-9, 'index_medium', 1.0);
beam.Nmax = 6;

x = linspace(-2.5, 2.5, 80)*1e-6;
y = x;
offset = [1.5e-6; 0; 0];

subplot(2, 3, 4);
data = beam.visualise('axis', 'y', 'range', {x, y});
imagesc(x, y, (data)), axis image;
drawNmax([0,0], beam.Nmax, 'r', ':');
drawNmax(-offset([2,1]).', beam.Nmax, 'w', '-');
caxr = caxis();
title('d');

subplot(2, 3, 5);
tbeam = beam.translateXyz(offset);
data = tbeam.visualise('axis', 'y', 'range', {x, y});
imagesc(x, y, (data)), axis image;
drawNmax(offset([2,1]).', beam.Nmax, 'r', ':');
drawNmax([0,0], tbeam.Nmax, 'w', '-');
caxis(caxr);
title('e');

subplot(2, 3, 6);
tbeam = beam.translateXyz(offset, 'Nmax', 20);
data = tbeam.visualise('axis', 'y', 'range', {x, y});
imagesc(x, y, (data)), axis image;
drawNmax(offset([2,1]).', beam.Nmax, 'r', ':');
drawNmax([0,0], tbeam.Nmax, 'k', '--');
caxis(caxr);
title('f');

% figure();
% obeam = tbeam.translateXyz(-offset, 'Nmax', beam.Nmax);
% data = beam.visualise('axis', 'y', 'range', {x, y});
% imagesc(x, y, (data)), axis image;


function drawNmax(cnt, nmax, color, type)

  ka = ott.utils.nmax2ka(nmax);
  a = ka ./ (2*pi/1064e-9);
  
  hold on;
  h = viscircles(cnt, a);
  set(h.Children(1), 'Color', color, 'LineStyle', type);
  delete(h.Children(2));
  hold off;
end