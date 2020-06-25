classdef Generic
% Constant definitions of different generic optical mediums.
%
% These mediums may not be physically accurate, nor are they accurate
% for all optical frequencies.  They are provided as a convenient set of
% defaults for modelling many common optical scenarios.
%
% Dielectric materials (pure real)
%   - Water       -- Non-absorbing water (index = 1.33)
%   - Polystyrene -- Approximation for polystyrene (index = 1.59)
%   - Glass       -- Approximation for some types of glass (index = 1.9)
%
% Conductive materials
%   - Gold        -- Conductive material for gold (index = 0.3 + 6.6i)
%
% Birefringent materials
%   - Vaterite    -- Birefringent material (index = [1.55, 1.55, 1.65])

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties (Constant)

    % Dielectric (pure real)
    Water = ott.beam.medium.Dielectric(1.33);
    Polystyrene = ott.beam.medium.Dielectric(1.59);
    Glass = ott.beam.medium.Dielectric(1.9);

    % Conductive
    Gold = ott.beam.medium.Dielectric(0.314 + 6.62i);

    % Birefringent
    Vaterite = ott.beam.medium.Dielectric(diag([1.55, 1.55, 1.65]));

  end

end

