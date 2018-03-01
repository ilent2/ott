% Example of the trapping_landscape code which could be used to produce 
% figure 3 in Nieminen et al., 2007, Journal of Optics. The nuts and bolts 
% of the how is included in trapping_landscape.m, this code simply scripts 
% for the parameters used in the article. More analysis of Trapping 
% landscapes and the properties of trapped microspheres appears in Stilgoe 
% et al., 2008, Optics Express. 
%
% The high resolution version takes about ~5000 seconds on a Core2 Duo 6600
% with 6GB RAM. The low resolution version takes about ~1030 seconds on the
% same machine.
%
% PACKAGE INFO

% Low res version.
size_range_rad=linspace(1e-2,3.25,100)*1e-6/2; %radius in SI
index_range=linspace(1.34,2.66,50);            %absolute refractive index

% % High res version.
% size_range_rad=linspace(1e-2,3.25,300)*1e-6/2; %radius in SI
% index_range=linspace(1.34,2.66,100);           %absolute refractive index

tic
%to see the proceedure for calculating a range of particle sizes and
%refractive indexes open trapping_landscape
structurelandscape=trapping_landscape(index_range,size_range_rad);
toc


%plot the minimum force. This is essentially the figure which appears in
%the article as figure 3.
tempmin=structurelandscape.minforce;
tempmin(tempmin>0)=NaN;
figure(1)
[d,pax]=contourf(size_range_rad*2*1e6,index_range,tempmin,20);
set(pax,'edgecolor','none')
xlabel('particle diameter [{{\mu}{m}}]');
ylabel('refractive index [unitless]');
grid on
cax=colorbar;
ylabel(cax,'Q_z^{min} [unitless]');

%plot the maximum force.
figure(2)
[d,pax]=contourf(size_range_rad*2*1e6,index_range,structurelandscape.maxforce,20);
set(pax,'edgecolor','none')
xlabel('particle diameter [{{\mu}{m}}]');
ylabel('refractive index [unitless]');
grid on
cax=colorbar;
ylabel(cax,'Q_z^{max} [unitless]');

%plot the z quilibrium position.
figure(3)
[d,pax]=contourf(size_range_rad*2*1e6,index_range,structurelandscape.zequilibrium,20);
set(pax,'edgecolor','none')
xlabel('particle diameter [{{\mu}{m}}]');
ylabel('refractive index [unitless]');
grid on
cax=colorbar;
ylabel(cax,'z [k{\lambda}]');

%plot the axial stiffness.
figure(4)
[d,pax]=contourf(size_range_rad*2*1e6,index_range,structurelandscape.zstiffness,20);
set(pax,'edgecolor','none')
xlabel('particle diameter [{{\mu}{m}}]');
ylabel('refractive index [unitless]');
grid on
cax=colorbar;
ylabel(cax,'k_z [k^{-1}\lambda^{-1}]');

%plot the transverse stiffness.
figure(5)
[d,pax]=contourf(size_range_rad*2*1e6,index_range,structurelandscape.xstiffness,20);
set(pax,'edgecolor','none')
xlabel('particle diameter [{{\mu}{m}}]');
ylabel('refractive index [unitless]');
grid on
cax=colorbar;
ylabel(cax,'k_r [k^{-1}\lambda^{-1}]');