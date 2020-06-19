classdef StokesLambPm < ott.drag.StokesStarShaped
% Calculate drag coefficients using Lamb series and point matching.
% Inherits from :class:`StokesStarShaped`
%
% Properties
%   - inverse     -- Calculated from `forward`
%   - forward     -- Drag tensor calculated using point matching.
%   - viscosity   -- Viscosity of medium
%   - shape       -- A star shaped particle describing the geometry
%
% Additional methods/properties inherited from :class:`Stokes`.
%
% Based on and incorporating code by Lachlan Gibson (2016)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    seriesOrder = 6
  end

  properties (SetAccess=protected, Hidden)
    rtp     % Particle shape (depends on shape)
    NM      % Coefficient list (depends on series order)
    A       % Coefficient matrix (depends on shape/series order)
  end

  properties (Dependent)
    forwardInternal       % Calculate forward drag tensor
  end

  methods
    function drag = StokesLambPm(varargin)
      % Construct star-shaped drag method using point-matching method.
      %
      % Usage
      %   drag = StokesLambPm(shape, ...)
      %
      % Parameters
      %   - shape (ott.shapes.Shape) -- Star shaped particle.
      %
      % Parameters can also be passed as named arguments.
      % Additional parameters passed to base.

      % Only need constructor for doc/help
      drag = drag@ott.drag.StokesStarShaped(varargin{:});
    end
  end

  methods (Hidden)
    function v = linearFlow(drag, xyz, direction)
      % Construct a linear flow field

      assert(numel(direction) == 3, 'direction must be 3 element vector');
      assert(size(xyz, 1) == 3 && ismatrix(xyz), 'xyz must be 3xN matrix');

      v = repmat(direction(:), 1, size(xyz, 2));
    end

    function v = rotationalFlow(drag, xyz, direction)
      % Construct a rotational flow field

      assert(numel(direction) == 3, 'direction must be 3 element vector');
      assert(size(xyz, 1) == 3 && ismatrix(xyz), 'xyz must be 3xN matrix');

      % 1x1xN
      X = reshape(xyz(1, :), 1, 1, []);
      Y = reshape(xyz(2, :), 1, 1, []);
      Z = reshape(xyz(3, :), 1, 1, []);
      Zero = zeros(size(X));

      % 3x3xN
      crossmat = [Zero, Z, -Y; -Z, Zero, X; Y, -X, Zero];

      % 3xN
      v = squeeze(sum(crossmat .* direction(:).', 2));
    end

    function rtp = surfacePoints(drag)
      % Calculate surface points for shape

      k = 20;
      n = drag.seriesOrder;

      % Generate grid of theta/phi
      N=[ceil(k*n/2+1),k*n]; % Angular grid size
      th=pi*linspace(0+0.5/N(1),1-0.5/N(1),N(1)); % Polar angles
      ph=2*pi*linspace(0,1-1/N(2),N(2)); % Azimuthal angles
      [PH, TH] = meshgrid(ph, th);

      % Calculate radial coordinates
      R = drag.shape.starRadii(TH(:), PH(:));

      % Package outputs
      rtp = [R(:), TH(:), PH(:)].';
    end

    function NM = coeffList(drag)
      % Generate coefficient list

      n = [0, drag.seriesOrder];

      % Number of coefficients (unknowns)
      n_co=3*n(1)*(n(1)+2)+3*n(2)*(n(2)+2);

      % Preallocate memory
      NM = zeros(5, n_co);

      ii=1;
      for ss=[1,2]
        for nn=1:n(ss)
          for mm=0:nn
            for tt=[1,2]
              for ll=1:3
                if ~(mm==0&&tt==2)
                  NM(1, ii)=nn;
                  NM(2, ii)=mm;
                  NM(3, ii)=ll;
                  NM(4, ii)=tt;
                  NM(5, ii)=ss;
                  ii=ii+1;
                end
              end
            end
          end
        end
      end

      NM = NM.';
    end

    function A = coeffMatrix(drag)
      % Evaluate StokesSphCoeff with this particle shape

      % Evaluate Associated Legendre values
      ALF = cell(1, drag.seriesOrder+1);
      for ii = 1:length(ALF)
        ALF{ii} = legendre(ii, cos(drag.rtp(2, :)), 'norm').';
      end

      Nrtp = size(drag.rtp, 2);
      rtp = drag.rtp;
      NM = drag.NM;
      rmax = max(drag.rtp(1, :));

      % Build coefficient matrix
      A = zeros(3*Nrtp, size(drag.NM, 1));
      for comp = 1:3

        % Calculate coefficient values (order must match solveSystem)
        Ac = StokesSphCoeffs(rtp(1, :).', rtp(2, :).', rtp(3, :).', ALF, ...
              drag.viscosity, comp, rmax);

        % Assemble coefficients
        for jj = 1:size(NM, 1)
          A((1:Nrtp) + (comp-1)*Nrtp, jj) = ...
            Ac{NM(jj, 3), NM(jj, 4), NM(jj, 5)}(NM(jj, 1),NM(jj, 2));
        end
      end
    end

    function sol = solveSystem(drag, vel)
      % Solve system of equations for velocity field

      assert(ismatrix(vel) && size(vel, 1) == 3, 'vel must be 3xN');

      % Construct RHS vector (transforming x and y velocities)
      B = [vel(1, :).*sin(drag.rtp(2, :)), ...
           vel(2, :).*sin(drag.rtp(2, :)), ...
           vel(3, :)].';

      % Solve system (A should be pre-computed)
      sol = drag.A\B; % Solve system

      % Remove Legendre normalisation from solution
      NM = drag.NM;
      F = (-1).^NM(:, 2).*sqrt((NM(:, 1)+0.5) ...
          .*factorial(NM(:, 1)-NM(:, 2))./factorial(NM(:, 1)+NM(:, 2)));
      sol = sol .* F;

      % Remove radial normalisation
      ind = NM(:, 5) == 1;
      rmax = max(drag.rtp(1, :));
      sol(ind) = sol(ind).*rmax.^-NM(ind, 1); % Remove radial normalisation
    end

    function col = calculateColumn(drag, velocities)
      % Calculate column of drag tensor and package results

      sol = drag.solveSystem(velocities);
      NM = drag.NM;

      % Force
      force=-4*pi*...
        [sol(NM(:,1)==1&NM(:,2)==1&NM(:,3)==1&NM(:,4)==1&NM(:,5)==2),...
        sol(NM(:,1)==1&NM(:,2)==1&NM(:,3)==1&NM(:,4)==2&NM(:,5)==2),...
        -sol(NM(:,1)==1&NM(:,2)==0&NM(:,3)==1&NM(:,4)==1&NM(:,5)==2)];

      % Torque
      torque=-8*pi*drag.viscosity*...
        [sol(NM(:,1)==1&NM(:,2)==1&NM(:,3)==3&NM(:,4)==1&NM(:,5)==2),...
        sol(NM(:,1)==1&NM(:,2)==1&NM(:,3)==3&NM(:,4)==2&NM(:,5)==2),...
        -sol(NM(:,1)==1&NM(:,2)==0&NM(:,3)==3&NM(:,4)==1&NM(:,5)==2)];

      col = [force(:); torque(:)];
    end

    function drag = setupProblem(drag)
      % Pre-compute coefficient list and matrix.
      %
      % This method is called by forward when rtp is not set.
      % Once this method is called, shape and seriesOrder are ignored
      % for the returned instance (original is unchanged).
      %
      % Usage
      %   drag = drag.setupProblem

      drag.rtp = drag.surfacePoints();
      drag.NM = drag.coeffList();
      drag.A = drag.coeffMatrix();  % Needs rtp and NM
    end
  end

  methods % Getters/setters
    function D = get.forwardInternal(drag)

      % Setup problem (if not already done)
      if isempty(drag.rtp)
        drag = drag.setupProblem();
      end

      % Calculate shape Cartesian coordinates
      xyz = ott.utils.rtp2xyz(drag.rtp);

      % Calculate drag tensor
      D = zeros(6, 6);
      D(:, 1) = drag.calculateColumn(drag.linearFlow(xyz, [1;0;0]));
      D(:, 2) = drag.calculateColumn(drag.linearFlow(xyz, [0;1;0]));
      D(:, 3) = drag.calculateColumn(drag.linearFlow(xyz, [0;0;1]));
      D(:, 4) = drag.calculateColumn(drag.rotationalFlow(xyz, [1;0;0]));
      D(:, 5) = drag.calculateColumn(drag.rotationalFlow(xyz, [0;1;0]));
      D(:, 6) = drag.calculateColumn(drag.rotationalFlow(xyz, [0;0;1]));
    end

    function drag = set.seriesOrder(drag, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'seriesOrder must be positive numeric scalar');
      drag.seriesOrder = val;
    end
  end
end
