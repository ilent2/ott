classdef AppProducer < handle
% Base class for applications which produce a single output.
%
% Properties (uninitialised)
%   - VariableName -- Output variable name widget
%   - UpdateButton -- Update button and auto-update toggle
%   - ShowPreviewCheckBox -- (Optional) Toggle to control preview status
%   - Data -- Data produced by the process.
%
% Abstract methods
%   - generateData
%
% Callback methods
%   - updateVariableNameCb -- Called when output variable name changes
%   - updateParametersCb -- Called when input parameters change
%   - previewToggledCb -- Called when the preview checkbox is changed
%   - updateButtonCb -- Called by the update button
%
% Methods
%   - updatePreview -- Called when the preview should be updated
%   - update -- Called when data should be updated
%   - writeData -- Called when data should be written

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (GetAccess=public, SetAccess=private)
    Data
  end

  % Widgets
  properties (Access=public)
    VariableName          ott.ui.support.OutputVariableEntry
    ShowPreviewCheckBox   matlab.ui.control.CheckBox
    UpdateButton          % ott.ui.support.UpdateButtonBase
  end
  
  % Callbacks
  methods (Access=private)
    function updateVariableNameCb(app, ~)
      % Called when output variable name changes
      app.writeData();
    end
    
    function previewToggledCb(app, ~)
      % Called when the preview checkbox is changed
      if app.ShowPreviewCheckBox.Value
        app.updatePreview();
      end
    end
    
    function updateButtonCb(app, ~)
      % Called by the update button
      app.update();
    end
  end
  
  methods (Access=protected)
    function updateParametersCb(app, ~)
      % Called when input parameters change
      app.update();
    end
  end
  
  methods (Abstract, Access=protected)
    generateData(app)
  end
  
  % Regular methods
  methods (Access=protected)
    function updatePreview(~)
      % Overload this methods to update the preview
    end
    
    function writeData(app)
      % Write data if available
      if ~isempty(app.VariableName.Value)
        app.VariableName.writeVariable(app.Data);
      end
    end
    
    function update(app)
      % Update the data (and output/preview)
      
      % Generate new data
      app.Data = app.generateData();
      
      % Write new data
      app.writeData();
      
      % Call preview method (if required)
      if ~isempty(app.Data) && ...
          (isempty(app.ShowPreviewCheckBox) || app.ShowPreviewCheckBox.Value)
        app.updatePreview();
      end
    end
  end
  
  methods
    function app = AppProducer()
      
      % Check required things are created
      assert(~isempty(app.VariableName), ...
        'VariableName must be created before calling constructor.');
      assert(~isempty(app.UpdateButton), ...
        'UpdateButton must be created before calling constructor.');
      
      % Connect callbacks
      app.VariableName.ValueChangedFcn = @(o,e) app.updateVariableNameCb(e);
      addlistener(app.UpdateButton, "UpdateCalled", ...
          @(h,e) app.updateButtonCb(e));
      if ~isempty(app.ShowPreviewCheckBox)
        app.ShowPreviewCheckBox.ValueChangedFcn = @(o,e) app.previewToggledCb(e);
      end
      
    end
  end
end