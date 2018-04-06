%Example farfield plotting code:

%build beam:
[a1,b1]=bsc_pointmatch_farfield(20,1,[ 1 1 lg_mode_w0([1,1],65) 1 1 -1i 90 0 0 0 ]);
[n,m]=combined_index(find(abs(a1)|abs(b1)));

%build grid:
nt=80;
[x,y,z]=sphere(nt);

%generate angular points for farfield:
[~,theta,phi]=xyz2rtp(x,y,z);

%find far-field in theta, phi:
[E,H]=farfield(n,m,a1,b1,0*a1,0*b1,theta(:),phi(:));


%% find radiant flux:
I=reshape(sum(abs(E).^2,2),[nt+1,nt+1]);

%% find circular polarisation: 
% horizontal and vertical require a coordinate transformation:
% (theta,phi)->(x,y)
R=reshape(abs(E(:,2)+1i*E(:,3)).^2,[nt+1,nt+1]); %can't remember which is left and right
L=reshape(abs(E(:,2)-1i*E(:,3)).^2,[nt+1,nt+1]);

%%
figure(1)
surf(x,y,z,R,'facecolor','interp','edgecolor','none')
axis equal
view([0,270])
title('right')
figure(2)
surf(x,y,z,L,'facecolor','interp','edgecolor','none')
axis equal
view([0,270])
title('left')