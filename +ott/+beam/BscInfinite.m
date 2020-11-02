classdef BscInfinite < ott.beam.BscBeam
% A Beam with stored Bsc data which grows as it is used.
% Inherits from :class:`ott.beam.BscBeam`.
%
% Useful for plane wave beams and other infinite extent beams.
%
% Properties
%   - data            -- Internal BSC instance describing beam
%
% Methods
%   - getData         -- Get data for specific Nmax
%   - recalculate     -- Update the internal data for new Nmax
%
% Abstract methods
%   - recalculateInternal     -- Update the internal data for new Nmax

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function varargout = getData(beam, Nmax)
      % Get BSC data for specific Nmax
      %
      % This function applies translations/rotations and grows the
      % beam data if it is not currently at least ``Nmax`` in size.
      %
      % Usage
      %   [data, beam] = beam.getData(Nmax)

      % Grow beam if required
      if beam.data.Nmax < Nmax
        beam = beam.recalculate(Nmax);
      end

      [varargout{1:nargout}] = getData@ott.beam.BscBeam(beam, Nmax);
    end

    function [beam, vswfData] = recalculate(beam, Nmax, varargin)
      % Re-calculate BSC data for specified Nmax.
      %
      % Only re-calculates the beam if Nmax exceeds current Nmax.
      %
      % Usage
      %   [beam, vswfData] = beam.recalculate(Nmax, ...)
      %
      % Optional named arguments
      %   - vswfData (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an emtpy VswfData structure.

      ott.utils.nargoutCheck(beam, nargout);

      if ~isempty(beam.data) && Nmax < beam.data.Nmax
        return;   % Nothing to do
      end

      p = inputParser;
      p.addParameter('vswfData', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.parse(varargin{:});

      % Defer to beam specific method
      [beam.data, vswfData] = beam.recalculateInternal(...
          Nmax, p.Results.vswfData);
    end
  end
end
