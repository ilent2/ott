% ott.optics.geometric Scattering calculations using geometric optics
%
% Geometric optics involves describing scattering particles by
% by surfaces and applying Snell's law to calculate reflected and
% refracted ray directions.
%
% Geometric optics uses the Ray beam, which is a specialisation of the
% PlaneWave beam with a finite area (and finite power).
% Scattering by a plane is equivalent to scattering of a plane wave by
% an infinite plane, except the scattered beams are :class:`Ray`s.
%
% For arbitrary shapes, scattering is implemented by calculating the
% intersection of rays with the shape and using the surface normals
% to calculate the scattering.
%
% In the future, specific shapes might get specialised implementations.
%
% Particles
%   Shape         -- Describes scattering by a ott.shapes.Shape instance
%
% Lenses
%   ThinLens      -- Thin lens approximation
%   SphericalLens -- Spherical lens approximation
%
% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file
