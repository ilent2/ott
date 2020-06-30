function [E, H, data] = emField(krtp, type, nm, ab, varargin)
% EMFIELD calculates field from the vector spherical wave functions
%
% [E, H, data] = emField(krtp, type, nm, ab, ...) calculates the
% E and H field for unit-less spherical coordinates krtp.
% krtp is a Nx3 matrix of [radial, polar, azimuthal] coordinates.
% type must be 'incoming', 'regular', or 'outgoing'.
%
% Optional named parameters:
%     'saveData'  bool    saves data that can be used for repeated
%         calculations of the fields at these locations (default: nargout==3).
%     'data'      data    data to use from previous calculation
%     'calcE'     bool    calculate E field (default: true)
%     'calcH'     bool    calculate H field (default: true)
%
% If internal fields are calculated only the theta and phi components
% of E are continuous at the boundary. Conversely, only the kr component of
% D is continuous at the boundary.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

import ott.utils.*

ip = inputParser;
ip.addParameter('calcE', true);
ip.addParameter('calcH', true);
ip.addParameter('saveData', nargout == 3);
ip.addParameter('data', []);
ip.addParameter('verbose', false);
ip.parse(varargin{:});

verbose = ip.Results.verbose;

ott.warning('internal');

% Unpack n and m vectors
n=nm(1:size(nm,1)/2,1);
m=nm(size(nm,1)/2+1:size(nm,1),1);

% Unpack a and b vectors
p=ab(1:size(ab,1)/2,1);
q=ab(size(ab,1)/2+1:size(ab,1),1);

% Get unique rtp (faster calculation)
[r_new,~,indR]=unique(krtp(:, 1));
[theta_new,~,indTheta]=unique(krtp(:, 2));
[phi_new,~,indPhi]=unique(krtp(:, 3));

% r_new=r_new;
r_new(r_new==0)=1e-15;
%ab length can be bigger than cd or pq as it is the beam we start with.

%look for biggest and smallest elements in the matrix, kill elements
%less than tol of the max.

if verbose
    disp(['emfieldxyz behaviour:', num2str(behaviour)]);
end

E = zeros(size(krtp));
H = E;

un=unique(n);

%NumberOfLoops=length(calcvecels)
if verbose
    tic
end

% Allocate memory for output data
data = [];
if ip.Results.saveData
  data = zeros(numel(indTheta), 0);
end

% Start a counter for accessing the data
if ~isempty(ip.Results.data)
  dataCount = 0;
end

for nn = 1:max(un)
    if verbose
        disp(['emfieldxyz nn:', num2str(nn)]);
    end
    Nn = 1/sqrt(nn*(nn+1));
    vv=find(n==nn);

    if ~isempty(vv)

      kr=r_new(indR);

      if isempty(ip.Results.data)

        [Y,Ytheta,Yphi] = spharm(nn,m(vv),theta_new,zeros(size(theta_new)));

        switch type
          case 'incoming'
            [hn,dhn]=sbesselh2(nn,r_new);
            hn = hn ./ 2;
            dhn = dhn ./ 2;

          case 'outgoing'
            [hn,dhn]=sbesselh1(nn,r_new);
            hn = hn ./ 2;
            dhn = dhn ./ 2;

          case 'regular'
            [hn,dhn]=sbesselj(nn,r_new);

          otherwise
            error('Unknown beam type');
        end

        [M,PHI]=meshgrid(1i*m(vv),phi_new);

        expimphi=exp(M.*PHI);

        hnU=hn(indR);
        dhnU=dhn(indR);

        % Create full Y, Ytheta, Yphi, expimphi matrices (opt, R2018a)
        expimphif = expimphi(indPhi, :);
        YExpf = Y(indTheta, :).*expimphif;
        YthetaExpf = Ytheta(indTheta, :).*expimphif;
        YphiExpf = Yphi(indTheta, :).*expimphif;

        % Save the data if requested
        if ip.Results.saveData
          data(:, end+1) = hnU;
          data(:, end+1) = dhnU;
          data(:, end+(1:size(Ytheta, 2))) = YExpf;
          data(:, end+(1:size(Ytheta, 2))) = YthetaExpf;
          data(:, end+(1:size(Ytheta, 2))) = YphiExpf;
        end

      else

        % Load the data if present
        hnU = ip.Results.data(:, dataCount+1);
        dataCount = dataCount + 1;
        dhnU = ip.Results.data(:, dataCount+1);
        dataCount = dataCount + 1;
        YExpf = ip.Results.data(:, dataCount+(1:length(vv)));
        dataCount = dataCount + length(vv);
        YthetaExpf = ip.Results.data(:, dataCount+(1:length(vv)));
        dataCount = dataCount + length(vv);
        YphiExpf = ip.Results.data(:, dataCount+(1:length(vv)));
        dataCount = dataCount + length(vv);

      end

      pidx = full(p(vv));
      qidx = full(q(vv));

      % Now we use full matrices, we can use matmul (opt, R2018a)
      if ip.Results.calcE
        E(:,1)=E(:,1)+Nn*nn*(nn+1)./kr.*hnU.*YExpf*qidx(:);
        E(:,2)=E(:,2)+Nn*(hnU.*YphiExpf*pidx(:) + dhnU.*YthetaExpf*qidx(:));
        E(:,3)=E(:,3)+Nn*(-hnU.*YthetaExpf*pidx(:) + dhnU.*YphiExpf*qidx(:));
      end

      if ip.Results.calcH
        H(:,1)=H(:,1)+Nn*nn*(nn+1)./kr.*hnU.*YExpf*pidx(:);
        H(:,2)=H(:,2)+Nn*((hnU(:).*YphiExpf)*qidx(:)+(dhnU(:).*YthetaExpf)*pidx(:));
        H(:,3)=H(:,3)+Nn*((-hnU(:).*YthetaExpf)*qidx(:)+(dhnU(:).*YphiExpf)*pidx(:));
      end
    end
end

H=-1i*H; %LOOK HERE TO FIX STUFF

ott.warning('external');
