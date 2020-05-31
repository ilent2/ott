classdef Paraxial < ott.beam.abstract.CastBoth ...
    & ott.beam.properties.VariablePower
% Common stem for Gaussian type beams.
% Inherits from :class:`CastBoth` and
% :class:`ott.beam.properties.VariablePower`.
%
% Abstract methods
%   - isGaussian      -- Returns true if the beam is a Gaussian
%
% Casts
%   - Beam                -- (Inherited) Uses vswf.Bsc/Coherent
%   - GaussianDavis5      -- Cast for Gaussian, uses abstract.Gaussian
%   - paraxial.Paraxial   -- Raises an error, must be implemented in sub
%   - paraxial.Gaussian         -- Cast for Gaussian, uses abstract.Gaussian
%   - paraxial.HermiteGaussian  -- Cast for Gaussian, uses abstract.Gaussian
%   - paraxial.LaguerreGaussian -- Cast for Gaussian, uses abstract.Gaussian
%   - vswf.Bsc              -- Raises an error, must be implemented in sub
%   - vswf.Gaussian         -- Cast for Gaussian, uses abstract.Gaussian
%   - vswf.LaguerreGaussian -- Cast for Gaussian, uses abstract.Gaussian
%   - vswf.HermiteGaussian  -- Cast for Gaussian, uses abstract.Gaussian
%   - vswf.InceGaussian     -- Cast for Gaussian, uses abstract.Gaussian
%   - abstract.Gaussian     -- Creates a Gaussian beam like the current beam
%   - abstract.LaguerreGaussian -- Cast for Gaussian, uses abstract.Gaussian
%   - abstract.HermiteGaussian  -- Cast for Gaussian, uses abstract.Gaussian
%   - abstract.InceGaussian     -- Cast for Gaussian, uses abstract.Gaussian
%
% Additional casts inherited from base.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.utils.VariablePower.likeProperties(other, args);
    end
  end

  methods (Hidden, Abstract)
    isGaussian(obj)
  end

  methods
    function beam = ott.beam.paraxial.Paraxial(varargin)
      error('Cast to paraxial not supported for this beam type');
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      error('Cast to Bsc not supported for this beam type');
    end

    function beam = ott.beam.abstract.Gaussian(beam, varargin)
      % Create a Gaussian beam like the current beam

      assert(isa(beam, 'ott.beam.abstract.Paraxial'), ...
          'First argument must be abstract.Paraxial beam');

      % Check if the beam isa Gaussian
      assert(beam.isGaussian(), 'Beam must be a Gaussian beam');

      beam = castHelper(@ott.beam.abstract.Gaussian.like, ...
          beam, varargin{:});
    end

    function beam = ott.beam.GaussianDavis5(beam, varargin)
      % Cast using abstract.Gaussian
      beam = castGaussian(beam);
      beam = castHelper(@ott.beam.GaussianDavis5.like, beam, varargin{:});
    end

    %
    % Paraxial casts
    %

    function beam = ott.beam.paraxial.Gaussian(beam, varargin)
      % Cast using abstract.Gaussian
      beam = castGaussian(beam);
      beam = castHelper(@ott.beam.paraxial.Gaussian.like, beam, varargin{:});
    end

    function beam = ott.beam.paraxial.HermiteGaussian(beam, varargin)
      % Cast using abstract.Gaussian
      beam = castGaussian(beam);
      beam = castHelper(@ott.beam.paraxial.HermiteGaussian.like, ...
          beam, 'mmode', 0, 'nmode', 0, varargin{:});
    end

    function beam = ott.beam.paraxial.LaguerreGaussian(beam, varargin)
      % Cast using abstract.Gaussian
      beam = castGaussian(beam);
      beam = castHelper(@ott.beam.paraxial.LaguerreGaussian.like, ...
          beam, 'lmode', 0, 'pmode', 0, varargin{:});
    end

    %
    % VSWF casts
    %

    function beam = ott.beam.vswf.Gaussian(beam, varargin)
      % Cast using abstract.Gaussian
      beam = castGaussian(beam);
      beam = castHelper(@ott.beam.vswf.Gaussian.like, ...
          beam, varargin{:});
    end

    function beam = ott.beam.vswf.LaguerreGaussian(beam, varargin)
      % Cast using abstract.Gaussian
      beam = castGaussian(beam);
      beam = castHelper(@ott.beam.vswf.LaguerreGaussian.like, ...
          beam, 'lmode', 0, 'pmode', 0, varargin{:});
    end

    function beam = ott.beam.vswf.HermiteGaussian(beam, varargin)
      % Cast using abstract.Gaussian
      beam = castGaussian(beam);
      beam = castHelper(@ott.beam.vswf.HermiteGaussian.like, ...
          beam, 'mmode', 0, 'nmode', 0, varargin{:});
    end

    function beam = ott.beam.vswf.InceGaussian(beam, varargin)
      % Cast using abstract.Gaussian
      error('Not yet implemented');
      beam = castGaussian(beam);
      beam = castHelper(@ott.beam.vswf.InceGaussian.like, ...
          beam, varargin{:});
    end

    %
    % Abstract casts
    %

    function beam = ott.beam.abstract.LaguerreGaussian(beam, varargin)
      % Cast using abstract.Gaussian
      beam = castGaussian(beam);
      beam = castHelper(@ott.beam.abstract.LaguerreGaussian.like, ...
          beam, 'lmode', 0, 'pmode', 0, varargin{:});
    end

    function beam = ott.beam.abstract.HermiteGaussian(beam, varargin)
      % Cast using abstract.Gaussian
      beam = castGaussian(beam);
      beam = castHelper(@ott.beam.abstract.HermiteGaussian.like, ...
          beam, 'mmode', 0, 'nmode', 0, varargin{:});
    end

    function beam = ott.beam.abstract.InceGaussian(beam, varargin)
      % Cast using abstract.Gaussian
      error('Not yet implemented');
      beam = castGaussian(beam);
      beam = castHelper(@ott.beam.abstract.InceGaussian.like, ...
          beam, varargin{:});
    end
  end

  methods (Access=protected)
    function beam = castGaussian(beam)
      % Cast to abstract.Gaussian helper

      assert(isa(beam, 'ott.beam.abstract.Beam'), ...
          'First argument must be a abstract.Beam');
      ott.utils.nargoutCheck(beam, nargout);

      % Check we have work to do
      if ~isa(beam, 'ott.beam.abstract.Gaussian')
        beam = ott.beam.abstract.Gaussian(beam);
      end
    end
  end
end
