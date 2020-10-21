classdef Pointmatch < ott.bsc.Bsc
% Provides Bsc constructors for point-matching from near/far fields.
% Inherits from :class:`Bsc`.
%
% Static methods
%   - FromNearfield   -- Construct from near-field points
%   - FromFarfield    -- Construct from far-field points

  methods (Static)
    function [beam, data] = FromNearfield(xyz, Exyz, ci, varargin)
      % Construct a beam using near-field point matching
      %
      % Usage
      %   [beam, data] = ott.bsc.Pointmatch.FromNearfield(xyz, Exyz, ci, ...)
      %
      % Parameters
      %   - xyz (3xN numeric) -- Locations for point matching
      %   - Exyz (3xN numeric) -- Field values for point matching
      %   - ci (numeric) -- Combed index Modes to include in point matching.
      %
      % Optional named parameters
      %   - basis (enum) -- Near-field basis, can be any of
      %     'incoming', 'regular', or 'outgoing'.  Default: 'regular'.
      %
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      p = inputParser;
      p.addParameter('basis', 'regular');
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.parse(varargin{:});

      % Generate basis set of beams
      Z = sparse(max(ci), numel(ci));
      I = sparse(ci, 1:numel(ci), ones(size(ci)), max(ci), numel(ci));
      vswfBasis = ott.bsc.Bsc([I, Z], [Z, I]);

      % Calculate coefficient matrix
      % Note: This doesn't have any assumptions about TEM fields
      [ourE, data] = vswfBasis.efield(xyz, 'data', p.Results.data, ...
          'basis', p.Results.basis);

      % Do point-matching step
      fab = reshape(ourE.vxyz, 3*size(xyz, 2), 2*numel(ci)) \ Exyz(:);

      % Package output
      beam = ott.bsc.Bsc(fab(1:end/2), fab(end/2+1:end));
    end

    function [beam, data] = FromFarfield(rtp, Ertp, ci, varargin)
      % Construct a beam using far-field point matching
      %
      % Usage
      %   [beam, data] = ott.bsc.Pointmatch.FromFarfield(rtp, Ertp, ci, ...)
      %
      % Parameters
      %   - rtp (2xN | 3xN numeric) -- Locations for point matching
      %
      %   - Ertp (2xN | 3xN numeric) -- Field values for point matching
      %     Ignores radial component if 3xN numeric input.
      %
      %   - ci (numeric) -- Combed index Modes to include in point matching.
      %
      % Optional named parameters
      %   - basis (enum) -- Near-field basis, can be any of
      %     'incoming', or 'outgoing'.  Default: 'incoming'.
      %
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      p = inputParser;
      p.addParameter('basis', 'incoming');
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.parse(varargin{:});

      % Generate basis set of beams
      Z = sparse(max(ci), numel(ci));
      I = sparse(ci, 1:numel(ci), ones(size(ci)), max(ci), numel(ci));
      vswfBasis = ott.bsc.Bsc([I, Z], [Z, I]);

      % Calculate coefficient matrix
      % Note: This doesn't have any assumptions about TEM fields
      [ourE, data] = vswfBasis.efarfield(rtp, 'data', p.Results.data, ...
          'basis', p.Results.basis);

      % Remove radial component
      if size(Ertp, 1) == 3
        Ertp = Ertp(2:3, :);
      end

      % Do point-matching step
      fab = reshape(ourE.vrtp(2:3, :), 2*size(rtp, 2), 2*numel(ci)) \ Ertp(:);

      % Package output
      beam = ott.bsc.Bsc(fab(1:end/2), fab(end/2+1:end));
    end
  end
end

