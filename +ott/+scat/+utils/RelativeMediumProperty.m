classdef RelativeMediumProperty
% Declares a relativeMedium property
%
% Does not implement a class constructor.  The `relativeMedium` property
% is not initialised and must be set in the sub-class constructor.
%
% Properties
%   - relativeMedium -- The relative medium property.
%
% Hidden methods
%   - validateMedium  -- Method called to validate the relativeMedium

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    relativeMedium % The encapsulated ott.beam.medium.Relative object
  end

  methods (Hidden)
    function val = validateMediu(~, val)
      % Method called to validate the medium.
      %
      % If your scattering method requires a specific material,
      % overload this method with the additional conditions or casts.

      assert(isa(val, 'ott.beam.medium.Relative'), ...
          'relativeMedium must be of type ''ott.beam.medium.Relative''');
    end
  end

  methods % Getters/setters
    function particle = set.relativeMedium(particle, val)
      particle.relativeMedium = particle.validateMedium(val);
    end
  end
end

