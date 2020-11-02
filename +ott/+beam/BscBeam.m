classdef BscBeam < ott.beam.Beam
% Base class for BSC beam classes.
% Inherits from :class:`+ott.+beam.Beam`.
%
% Changing any of the beam data related properties causes the beam data
% to be cleared.  The next call to :meth:`getData` re-calculates the beam
% data.
%
% Properties
%   - data          -- Internal BSC instance describing beam
%
% Abstract methods
%   - recalculateInternal -- Beam specific method called by ``recalculate``.
%
% Methods
%   - getData       -- Get data for specific Nmax
%   - recalculate   -- Recalculate the beam data

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    data        % Internal BSC instance describing beam
  end

  methods (Abstract, Hidden)
    recalculateInternal(~)
  end

  methods
    function [beam, vswfData] = recalculate(beam, ~, varargin)
      % Re-calculate BSC data
      %
      % Only re-calculate the beam if there is no valid BSC data.
      %
      % Usage
      %   [beam, vswfData] = beam.recalculate(Nmax, ...)
      %
      % Parameters
      %   - Nmax -- Ignored.
      %
      %   - vswfData (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an emtpy VswfData structure.

      ott.utils.nargoutCheck(beam, nargout);

      p = inputParser;
      p.addParameter('vswfData', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.parse(varargin{:});

      if isempty(beam.data)
        [beam.data, vswfData] = beam.recalculateInternal(...
            [], p.Results.vswfData);
      end
    end

    function [data, beam] = getData(beam, Nmax)
      % Get BSC data for specific Nmax
      %
      % This function applies the translations/rotations for a
      % specific ``Nmax``.
      %
      % Usage
      %   [data, beam] = beam.getData(Nmax)

      % If no valid data, re-calculate data
      if isempty(beam.data)
        beam = beam.recalculate(Nmax);
      end

      % Apply rotation
      if all(all(beam.rotation ~= eye(3)))
        beam = beam.applyRotation('Nmax', Nmax);
      end

      % Apply translation
      if all(beam.translation ~= [0;0;0])
        beam = beam.applyTranslation('Nmax', Nmax);
      end

      data = beam.data;
    end
  end
end
