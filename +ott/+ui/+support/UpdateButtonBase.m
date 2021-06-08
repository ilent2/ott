classdef UpdateButtonBase < handle
% Base class for update button

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  events (NotifyAccess = protected)
    UpdateCalled    % Emitted when button clicked or auto-update enabled
  end
  
  properties (Abstract)
    AutoUpdate
  end
  
  properties (Dependent)
    UpdateCalledFcn   % Convinience method to set UpdateCalled callback
  end
  
  methods (Abstract)
    setError(app)
    setWarning(app)
    clearErrors(app)
  end
  
  methods
    function set.UpdateCalledFcn(obj, fcn)
      addlistener(obj, "UpdateCalled", fcn);
    end
  end
end