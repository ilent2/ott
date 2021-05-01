classdef (Abstract) AppBase < ott.ui.support.AppTwoColumn

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (SetAccess=protected)
    beam
  end
  
  properties (Access=public)
    PreviewUIAxes                         matlab.ui.control.UIAxes
  end
  
  methods (Access=protected)
    function updateBeamPreview(app)
      
      if isempty(app.beam)
        return;
      end
      
      app.beam.visNearfield('plot_axes', app.PreviewUIAxes, ...
        'axis', 'y', 'range', [1,1]*2e-6, 'field', 'Re(Ex)', ...
        'size', [60, 60]);
      
    end
    
    function createBeamPreview(app)

      % Create UIAxes
      app.PreviewUIAxes = uiaxes(app.RightPanel);
      title(app.PreviewUIAxes, 'Preview')
      xlabel(app.PreviewUIAxes, '')
      ylabel(app.PreviewUIAxes, '')
      app.PreviewUIAxes.XAxisLocation = 'origin';
      app.PreviewUIAxes.XTick = [];
      app.PreviewUIAxes.YAxisLocation = 'origin';
      app.PreviewUIAxes.YTick = [];
      app.PreviewUIAxes.Position = [10 10 373 app.windowSize(2)-20];
      
    end
  end
  
end