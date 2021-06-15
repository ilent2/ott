classdef FileSelector < ott.ui.support.LabeledWidget
% Creates a file entry with a browse button

% Copyright 2020 Isaac Lenton.
% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    EditField   matlab.ui.control.EditField
    Button      matlab.ui.control.Button
    Filter (:,:) cell
    Title (1, :) char
    PostSelect  matlab.ui.Figure
  end
  
  properties (Dependent)
    Value
    ValueChangedFcn
  end
  
  methods (Access=private)
    function BrowseClickedButtonCb(obj, ~)
      % Callback for browse button click
      
      % Display a file chooser dialog
      [file, path] = uigetfile(obj.Filter, obj.Title);

      % Put the GUI back on top
      if ishandle(obj.PostSelect)
        figure(obj.PostSelect);
      end

      % Check if the user selected a file
      if ~isequal(file,0)

          % Update the filename field and process the event
          obj.EditField.Value = fullfile(path, file);

          % Run the value changed callback
          obj.ValueChangedFcn([]);
      end
    end
  end

  methods
    function obj = FileSelector(parent, varargin)
      % Create a new file selector
      
      p = inputParser;
      p.addParameter('label', 'Spinner');
      p.addParameter('visible', 'on');
      p.parse(varargin{:});
      
      obj = obj@ott.ui.support.LabeledWidget(parent, ...
          'label', p.Results.label, 'visible', p.Results.visible);
        
      % Configure grid for additional entry
      obj.Grid.ColumnWidth = {'1x', 80, 20};
      
      % Create text field
      obj.EditField = uieditfield(obj.Grid, 'text');
      obj.EditField.Layout.Column = 2;
      obj.EditField.Layout.Row = 1;
      
      % Create button
      obj.Button = uibutton(obj.Grid, 'push');
      obj.Button.Text = '...';
      obj.Button.Layout.Column = 3;
      obj.Button.Layout.Row = 1;
      obj.Button.ButtonPushedFcn = @(h,e) obj.BrowseClickedButtonCb(e);
      
    end
  end
  
  methods
    function val = get.Value(obj)
      val = obj.EditField.Value;
    end
    
    function set.Value(obj, val)
      obj.EditField.Value = val;
    end
    
    function set.ValueChangedFcn(obj, val)
      obj.EditField.ValueChangedFcn = val;
    end
    
    function val = get.ValueChangedFcn(obj)
      val = obj.EditField.ValueChangedFcn;
    end
  end
end