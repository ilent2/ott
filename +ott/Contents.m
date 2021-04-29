% OTT Optical Tweezers Toolbox
%
% The optical tweezers toolbox is a collection of functions and classes
% that can be used to simulate particles held in optical traps.
% The toolbox is split into several sub-packages, listed bellow.
%
% To use a class from any of these packages, simply prefix the class
% name with the package name(s), for example::
%
%   beam = ott.beam.Gaussian()
%
% Sub-packages
%   +beam       - High level descriptions of optical tweezers beams
%   +bsc        - Low level methods for working with beam shape coefficients
%   +drag       - Methods for calculating drag tensors
%   +dynamics   - Tools for simulating particle dynamics
%   +particle   - High level particle descriptions (drag/tmatrix/mass/geometry)
%   +shape      - Methods for describing particle geometry
%   +tmatrix    - T-matrix description of optical scattering
%   +tools      - Useful tools for simulating dynamics and finding traps
%   +utils      - Utility functions (mostly used by other toolbox functions)
%   +ui         - Graphical user interface
%
% Copyright 2018-2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

