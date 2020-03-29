% Look at the angular_scaling result for BscPmGauss

figure();
subplot(1, 2, 1);
beam = ott.BscPmGauss('angular_scaling', 'sintheta');
beam.basis = 'incoming';
beam.visualiseFarfield('dir', 'neg');

subplot(1, 2, 2);
beam = ott.BscPmGauss('angular_scaling', 'tantheta');
beam.basis = 'incoming';
beam.visualiseFarfield('dir', 'neg');