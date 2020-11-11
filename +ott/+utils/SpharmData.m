classdef SpharmData
% Stores pre-computed VSWF field data, useful for repeated field calculations.
%
% Properties
%   - uphi            -- Unique list of phi for this data set
%   - utheta          -- Unique list of theta for this data set
%   - uci             -- Combined index for stored values
%
%   - expimphi        -- Pre-computed exp(1i*m.*phi)
%
%   - calculated      -- Tracks which spherical harmonics have been computed
%   - Y               -- Pre-computed scalar spherical harmonic
%   - Ytheta          -- Pre-computed theta angular derivative
%   - Yphi            -- Pre-computed phi angular derivative
%
% Methods
%   - VswfData        -- Constructor
%   - evaluate        -- Evaluate Y, Ytheta, and Yphi
%   - updateCi        -- Update mode indices in data set
%   - updatePhi       -- Update phi values in data set
%   - updateTheta     -- Update theta values in data set

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    uphi             % Unique list of phi for this data set
    utheta           % Unique list of theta for this data set
    um               % m index for stored values (expimphi values)
    uci              % Combined index for stored values (Y values)

    expimphi         % Pre-computed exp(1i*m.*phi)
    calculated       % Tracks which spherical harmonics have been computed
    Y                % Pre-computed scalar spherical harmonic
    Ytheta           % Pre-computed theta angular derivative
    Yphi             % Pre-computed phi angular derivative
  end

  methods
    function data = SpharmData()
      % Initialise and empty data set

      data.uphi = [];
      data.utheta = [];
      data.uci = [];

      data.expimphi = [];
      data.calculated = [];
      data.Y = [];
      data.Ytheta = [];
      data.Yphi = [];
    end

    function [Y, Ytheta, Yphi, data] = evaluate(data, ci, theta, phi)
      % Evaluate scalar harmonics and update stored data
      %
      % Usage
      %   [Y, Ytheta, Yphi, data] = data.evaluateYtp(ci, theta, phi)
      %
      % Parameters
      %   - ci (N numeric) -- Combined indices for modes to evaluate.
      %   - theta, phi (M numeric) -- Angular coordinates for point.
      %
      % Returns
      %   - Y, Ytheta, Yphi (NxM numeric) -- Spherical harmonics/derivatives.

      assert(isvector(ci) && isnumeric(ci), ...
          'ci must be numeric vector');
      assert(isvector(theta) && isnumeric(theta), ...
          'theta must be numeric vector');
      assert(isvector(phi) && isnumeric(phi), ...
          'phi must be numeric vector');

      % Update internal data for theta/phi/ci
      data = data.updateCi(ci);
      data = data.updateTheta(theta);
      data = data.updatePhi(phi);
      
      [phi, theta] = ott.utils.matchsize(phi(:), theta(:));
      keep = ~isnan(phi) & ~isnan(theta);

      % Get indices for elements
      [~, dci] = ismember(ci, data.uci);
      [~, dtheta] = ismember(theta, data.utheta);
      [~, dphi] = ismember(phi, data.uphi);
      [~, m] = ott.utils.combined_index(ci);
      [~, dm] = ismember(m, data.um);

      % Update elements we actually need
      data = data.updateSpharm(dci, dtheta(keep));
      
      % Allocate memory for output
      Y = zeros(numel(ci), numel(phi));
      Ytheta = Y;
      Yphi = Y;
      
      % Assign nan to invalid values
      Y(:, ~keep) = nan;
      Ytheta(:, ~keep) = nan;
      Yphi(:, ~keep) = nan;

      % Retrieve stored data
      if ~isempty(Y)
        oexpimphi = data.expimphi(dm, dphi(keep));
        Y(:, keep) = data.Y(dci, dtheta(keep)) .* oexpimphi;
        Ytheta(:, keep) = data.Ytheta(dci, dtheta(keep)) .* oexpimphi;
        Yphi(:, keep) = data.Yphi(dci, dtheta(keep)) .* oexpimphi;
      end
    end

    function data = updateSpharm(data, dci, dtheta)
      % Update elements of spharm that are not already computed
      %
      % Usage
      %   data = data.updateSpharm(rowidx, colidx)

      ci = data.uci(dci);
      [nn, mm] = ott.utils.combined_index(ci);
      theta = data.utheta(dtheta);

      % Calculate block with similar N values.
      % This is a little redundant, some values might be recalculated
      % if we are using lots of scattered coordinates.
      for ii = unique(nn(:)).'

        idx = nn == ii;
        odci = dci(idx);
        omm = mm(idx);

        % Find values that have been computed
        isCalculated = data.calculated(odci, dtheta);

        % Calculate
        if ~all(all(isCalculated))
          
          % Find uncalculated values
          midx = ~all(isCalculated, 2);
          tidx = ~all(isCalculated, 1);

          % Calculate
          [oY, oYtheta, oYphi] = ott.utils.spharm(ii, omm(midx), ...
              theta(tidx), zeros(size(theta(tidx))));

          % Store new values
          data.Y(odci(midx), dtheta(tidx)) = oY.';
          data.Ytheta(odci(midx), dtheta(tidx)) = oYtheta.';
          data.Yphi(odci(midx), dtheta(tidx)) = oYphi.';
          data.calculated(odci(midx), dtheta(tidx)) = true;
        end

      end
    end

    function data = updateCi(data, ci)
      % Ensure specified ci are in data set, otherwise compute them
      %
      % Usage
      %   data = data.updateCi(ci)

      ott.utils.nargoutCheck(data, nargout);

      % Find which um and uci need to be updated
      ci = unique(ci);
      [~, dci] = ismember(ci, data.uci);
      ci_new = ci(dci == 0);
      [~, m] = ott.utils.combined_index(ci_new);
      m = unique(m);
      [~, dm] = ismember(m, data.um);
      m_new = m(dm == 0);

      % Update expimphi
      if ~isempty(m_new)
        data.um = [data.um, m_new(:).'];
        if ~isempty(data.uphi)
          data.expimphi = [data.expimphi; exp(1i*m_new(:) .* data.uphi)];
        end
      end

      % Update Y/dtY/dpY
      if ~isempty(ci_new)

        % Allocate space, do actual calculation later
        data.uci = [data.uci, ci_new(:).'];
        data.Y = [data.Y; zeros(numel(ci_new), size(data.Y, 2))];
        data.Ytheta = [data.Ytheta; zeros(numel(ci_new), size(data.Y, 2))];
        data.Yphi = [data.Yphi; zeros(numel(ci_new), size(data.Y, 2))];
        data.calculated = [data.calculated; ...
            false(numel(ci_new), size(data.Y, 2))];
      end
    end

    function data = updatePhi(data, phi)
      % Ensure specified phi are in data set, otherwise calculate them
      %
      % Usage
      %   data = data.updatePhi(phi)

      ott.utils.nargoutCheck(data, nargout);

      % Find phis to update
      phi(isnan(phi)) = [];
      phi = unique(phi);
      [~, dphi] = ismember(phi, data.uphi);
      new_phi = phi(dphi == 0);

      % Update expimphi
      if ~isempty(new_phi) && ~isempty(data.um)
        data.uphi = [data.uphi, new_phi(:).'];
        data.expimphi = [data.expimphi, exp(1i*data.um.' .* new_phi(:).')];
      end

    end

    function data = updateTheta(data, theta)
      % Ensure specified theta are in data set, otherwise calculate them
      %
      % Usage
      %   data = data.updateTheta(theta)

      ott.utils.nargoutCheck(data, nargout);

      % Find phis to update
      theta(isnan(theta)) = [];
      theta = unique(theta);
      [~, dtheta] = ismember(theta, data.utheta);
      new_theta = theta(dtheta == 0);

      % Update spharm
      if ~isempty(new_theta)

        % Allocate space, do actual calculation later
        data.utheta = [data.utheta, new_theta(:).'];
        data.Y = [data.Y, zeros(size(data.Y, 1), numel(new_theta))];
        data.Ytheta = [data.Ytheta, zeros(size(data.Y, 1), numel(new_theta))];
        data.Yphi = [data.Yphi, zeros(size(data.Y, 1), numel(new_theta))];
        data.calculated = [data.calculated, ...
            false(size(data.Y, 1), numel(new_theta))];
      end
    end
  end
end

