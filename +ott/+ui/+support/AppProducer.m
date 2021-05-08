classdef AppProducer < handle
% Base class for applications which produce a single output.
%
% Properties (uninitialised)
%   - VariableName -- Output variable name widget
%   - UpdateButton -- Update button and auto-update toggle
%   - PreviewCheckBox -- (Optional) Toggle to control preview status
%   - Data -- Data produced by the process.
%
% Abstract methods
%   - GenerateData
%
% Callback methods
%   - UpdateVariableNameCb -- Called when output variable name changes
%   - UpdateParametersCb -- Called when input parameters change
%   - PreviewToggledCb -- Called when the preview checkbox is changed
%   - UpdateButtonCb -- Called by the update button
%
% Methods
%   - UpdatePreview -- Called when the preview should be updated
%   - Update -- Called when data should be updated
%   - WriteData -- Called when data should be written

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
    UpdateButton          ott.ui.support.UpdateButtonBase
  end
  
  % Callbacks
  methods (Access=private)
    function UpdateVariableNameCb(app, ~)
      % Called when output variable name changes
      app.WriteData();
    end
    
    function PreviewToggledCb(app, ~)
      % Called when the preview checkbox is changed
      if app.ShowPreviewCheckBox.Value
        app.UpdatePreview();
      end
    end
    
    function UpdateButtonCb(app, ~)
      % Called by the update button
      app.Update();
    end
  end
  
  methods (Access=protected)
    function UpdateParametersCb(app, ~)
      % Called when input parameters change
      app.Update();
    end
  end
  
  methods (Abstract, Access=protected)
    GenerateData(app)
  end
  
  % Regular methods
  methods (Access=protected)
    function UpdatePreview(~)
      % Overload this methods to update the preview
    end
    
    function WriteData(app)
      % Write data if available
      if ~isempty(app.VariableName.Value)
        app.VariableName.WriteVariable(app.Data);
      end
    end
    
    function Update(app)
      % Update the data (and output/preview)
      
      % Generate new data
      app.Data = app.GenerateData();
      
      % Write new data
      app.WriteData();
      
      % Call preview method (if required)
      if isempty(app.ShowPreviewCheckBox) || app.ShowPreviewCheckBox.Value
        app.UpdatePreview();
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
      app.VariableName.ValueChangedFcn = createCallbackFcn(app, ...
          @UpdateVariableNameCb, true);
      addlistener(app.UpdateButton, "UpdateCalled", ...
          @(h,e) app.UpdateButtonCb(e));
      if ~isempty(app.ShowPreviewCheckBox)
        app.ShowPreviewCheckBox.ValueChangedFcn = createCallbackFcn(app, ...
            @PreviewToggledCb, true);
      end
      
    end
  end
end