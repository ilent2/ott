classdef Method
% Base class for optical tweezers dynamics methods.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    temperature         % Temperature for stochastic simulations
  end

  methods
    function x = dynamics(obj)

      % Dynamic step size?

      for ii = 1:numt
        x = x + obj.step();
      end

    end
  end

end
