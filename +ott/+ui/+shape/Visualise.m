classdef Visualise < ott.ui.support.AppTwoColumn ...
    & ott.ui.support.GenerateCodeMenu ...
    & ott.ui.support.RefreshInputsMenu
% Generate a visualisation of a OTT geometric shape.
%
% This GUI can be launched from the launcher under
% Shape -> Visualise or running the following command:
%
%   ott.ui.shape.Visualise()

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Shape as input argument

  properties (Constant)
    cnameText = 'Visualise';

    nameText = 'Visualise Shape';

    aboutText = ['Generate a visualisation of an existing shape.  ', ...
      'This GUI can be used to create a surface or voxel based visualisation', ...
      'of an OTT shape.'];
    
    helpText = {ott.ui.shape.Visualise.aboutText, ...
      '', ...
      ['This GUI requires that there is a ott.shape.Shape object already ', ...
      'in the current matlab workspace.  If you don''t yet have a shape, ', ...
      'you can create one with the Shape>Simple GUI.'], ...
      '', ...
      ['shape - Specify the variable name of the OTT shape to visualise. ', ...
      'If this list is empty, try clicking ''File > Refresh Inputs'' or ', ...
      'check that the variable is a valid ott.shape.Shape object.'], ...
      '', ...
      ['Visualisation - Specify the type of visualisation to create.  ', ...
      'See the documentation for visualisation specific parameters.'], ...
      '', ...
      };
    
    windowName = ott.ui.shape.Visualise.nameText;
    windowSize = [640, 350];
  end

  % Properties that correspond to app components
  properties (Access = public)
    
    % Left panel
    MainGrid                    matlab.ui.container.GridLayout
    ShapeDropDown               ott.ui.support.VariableDropDown
    VisualisationDropDown       ott.ui.support.LabeledDropDown
    OriginDropDown              ott.ui.support.LabeledDropDown
    UpdateButton                ott.ui.support.UpdateCheckButton
    
    % Surface widgets
    SurfNormals                 matlab.ui.control.CheckBox
    SurfNormalScale             ott.ui.support.LabeledSpinner
    SurfEdges                   matlab.ui.control.CheckBox
    SurfSmooth                  matlab.ui.control.CheckBox
    
    % Voxel widgets
    VoxSpacingSpinner           ott.ui.support.LabeledSpinner
    VoxScaleSpinner             ott.ui.support.LabeledSpinner
    VoxEvenRange                matlab.ui.control.CheckBox
    
    % Right panel
    LoadingText                 matlab.ui.control.Label
    UIAxes                      matlab.ui.control.UIAxes
  end

  % Properties that correspond to apps with auto-reflow
  properties (Access = private)
    onePanelWidth = 576;
  end

  % Component initialization
  methods (Access = protected)
    
    function setDefaultValues(app, ~)
      
      % Generics
      app.ShapeDropDown.Value = '';
      app.VisualisationDropDown.Value = 'Surface';
      app.UpdateButton.AutoUpdate = true;
      app.OriginDropDown.Value = 'local';
      
      % Vox widgets
      app.VoxScaleSpinner.Value = 1;
      app.VoxSpacingSpinner.Value = 0.1;
      app.VoxEvenRange.Value = false;
      
      % Surf widgets
      app.SurfNormals.Value = false;
      app.SurfNormalScale.Value = 0.1;
      app.SurfEdges.Value = true;
      app.SurfSmooth.Value = false;
      
      % Update visible widgets
      app.UpdateVisibleWidgets();
      
      % Clear figure
      cla(app.UIAxes);
      title(app.UIAxes, '');
    end
    
    function code = generateCode(app)
      % Generate visualisation code
      code = {};
      
      code{end+1} = 'figure()';
      
      % Update shape preview
      switch app.VisualisationDropDown.Value
        case 'Surface'
          code{end+1} = 's = shape.surf(...';
          code{end+1} = ['  ''normalScale'', ' num2str(app.SurfNormalScale.Value) ', ...'];
          code{end+1} = ['  ''showNormals'', ' num2str(app.SurfNormals.Value) ', ...'];
          code{end+1} = ['  ''origin'', ''' app.OriginDropDown.Value ''');'];
          
          if ~app.SurfEdges.Value
            code{end+1} = '';
            code{end+1} = '% Hide edges';
            code{end+1} = 's.EdgeColor = ''none'';';
          end
          
          if app.SurfSmooth.Value
            code{end+1} = '';
            code{end+1} = '% Smooth lighting';
            code{end+1} = 'camlight(app.UIAxes);';
            code{end+1} = 'lighting(app.UIAxes, ''gouraud'');';
          end
          
        case 'Voxels'
          
          code{end+1} = 'plotoptions = {...';
          code{end+1} = '  ''MarkerFaceColor'', ''w'', ...';
          code{end+1} = '  ''MarkerEdgeColor'', [.5 .5 .5], ...';
          code{end+1} = ['  ''MarkerSize'', ' num2str(app.VoxScaleSpinner.Value) '};'];
          code{end+1} = '';
          
          code{end+1} = 'shape.voxels(...';
          code{end+1} = ['  ''spacing'', ' num2str(app.VoxSpacingSpinner.Value) ', ...'];
          code{end+1} = '  ''plotoptions'', plotoptions, ...';
          code{end+1} = ['  ''even_range'', ' num2str(app.VoxEvenRange.Value) ', ...'];
          code{end+1} = ['  ''origin'', ''' app.OriginDropDown.Value ''');'];
          
        otherwise
          error('Internal error');
      end
      
      code{end+1} = '';
      
      code{end+1} = '% Axis Labels';
      code{end+1} = 'xlabel(''X'');';
      code{end+1} = 'ylabel(''Y'');';
      code{end+1} = 'zlabel(''Z'');';
    end
    
    function UpdateVisibleWidgets(app, ~)
      
      rheight = 32;
      rng = 4:10;
      
      % Hide all widgets
      app.MainGrid.RowHeight(rng) = repmat({0}, 1, length(rng));
      
      % Disable all widgets
      widgets = findobj(app.MainGrid.Children(rng), '-property', 'Enable');
      set(widgets, 'Enable', 'off');
      
      switch app.VisualisationDropDown.Value
        case 'Surface'
          app.MainGrid.RowHeight{app.SurfEdges.Layout.Row} = rheight;
          app.MainGrid.RowHeight{app.SurfNormals.Layout.Row} = rheight;
          app.MainGrid.RowHeight{app.SurfNormalScale.Layout.Row} = rheight;
          app.MainGrid.RowHeight{app.SurfSmooth.Layout.Row} = rheight;
          app.SurfEdges.Enable = 'on';
          app.SurfNormals.Enable = 'on';
          app.SurfNormalScale.Enable = app.SurfNormals.Value;
          app.SurfSmooth.Enable = 'on';
        case 'Voxels'
          app.MainGrid.RowHeight{app.VoxScaleSpinner.Layout.Row} = rheight;
          app.MainGrid.RowHeight{app.VoxSpacingSpinner.Layout.Row} = rheight;
          app.MainGrid.RowHeight{app.VoxEvenRange.Layout.Row} = rheight;
          app.VoxScaleSpinner.Enable = 'on';
          app.VoxSpacingSpinner.Enable = 'on';
          app.VoxEvenRange.Enable = 'on';
        otherwise
          error('Internal error');
      end
    end
    
    function ValueChangedCb(app, ~)
      if app.UpdateButton.AutoUpdate
        app.UpdateCb();
      end
    end
    
    function VisualiseChangedCb(app, ~)
      app.UpdateVisibleWidgets();
      app.ValueChangedCb();
    end
    
    function SurfNormalsValueChangedCb(app, ~)
      % Callback for surface normals toggle
      app.SurfNormalScale.Enable = app.SurfNormals.Value;
      app.ValueChangedCb();
    end
    
    function UpdateCb(app, ~)
      
      % Check for work to do
      if isempty(app.ShapeDropDown.Value)
        return;
      end
      
      % Get shape
      shape = app.ShapeDropDown.Variable;
      
      % Display loading text
      app.LoadingText.Visible = 'on';
      drawnow nocallbacks;
      
      % Update shape preview
      switch app.VisualisationDropDown.Value
        case 'Surface'
          s = shape.surf('axes', app.UIAxes, ...
            'normalScale', app.SurfNormalScale.Value, ...
            'showNormals', app.SurfNormals.Value, ...
            'origin', app.OriginDropDown.Value);
          
          % Hide edges
          if ~app.SurfEdges.Value
            s.EdgeColor = 'none';
          end
          
          % Smooth lighting
          if app.SurfSmooth.Value
            camlight(app.UIAxes);
            lighting(app.UIAxes, 'gouraud');
          end
          
          % No title for surf (at the moment)
          title(app.UIAxes, '');
          
        case 'Voxels'
          
          plotoptions = {...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', [.5 .5 .5], ...
            'MarkerSize', app.VoxScaleSpinner.Value};
          
          shape.voxels('axes', app.UIAxes, ...
            'spacing', app.VoxSpacingSpinner.Value, ...
            'plotoptions', plotoptions, ...
            'even_range', app.VoxEvenRange.Value, ...
            'origin', app.OriginDropDown.Value);
        otherwise
          error('Internal error');
      end
      
      % Labels
      xlabel(app.UIAxes, 'X');
      ylabel(app.UIAxes, 'Y');
      zlabel(app.UIAxes, 'Z');
      
      app.LoadingText.Visible = 'off';
      
    end

    % Create UIFigure and components
    function createLeftComponents(app)
      
      % Left layout grid
      app.MainGrid = uigridlayout(app.LeftPanel);
      app.MainGrid.RowHeight = repmat({32}, 1, 12);
      app.MainGrid.RowHeight{end-1} = '1x';
      app.MainGrid.ColumnWidth = {'1x'};
      app.MainGrid.RowSpacing = 0;

      % Create Shape variable select
      app.ShapeDropDown = ott.ui.support.VariableDropDown(app.MainGrid, ...
        'label', 'Shape', 'filter', 'ott.shape.Shape');
      app.ShapeDropDown.Layout.Row = 1;
      app.ShapeDropDown.Layout.Column = 1;
      app.registerRefreshInput(app.ShapeDropDown);
      app.ShapeDropDown.ValueChangedFcn = createCallbackFcn(app, ...
        @ValueChangedCb, true);
      
      % Create visualisation drop down
      app.VisualisationDropDown = ott.ui.support.LabeledDropDown(app.MainGrid, ...
        'label', 'Visualisation');
      app.VisualisationDropDown.Items = {'Surface', 'Voxels'};
      app.VisualisationDropDown.Layout.Row = 2;
      app.VisualisationDropDown.Layout.Column = 1;
      app.VisualisationDropDown.ValueChangedFcn = createCallbackFcn(app, ...
        @VisualiseChangedCb, true);
      
      % Create visualisation drop down
      app.OriginDropDown = ott.ui.support.LabeledDropDown(app.MainGrid, ...
        'label', 'Origin');
      app.OriginDropDown.Items = {'Local', 'Global'};
      app.OriginDropDown.ItemsData = {'local', 'global'};
      app.OriginDropDown.Layout.Row = 3;
      app.OriginDropDown.Layout.Column = 1;
      app.OriginDropDown.ValueChangedFcn = createCallbackFcn(app, ...
        @ValueChangedCb, true);
      
      % Vox: Spacing
      app.VoxSpacingSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
        'label', 'Voxel Spacing');
      app.VoxSpacingSpinner.Layout.Row = 4;
      app.VoxSpacingSpinner.Layout.Column = 1;
      app.VoxSpacingSpinner.Limits = [0, Inf];
      app.VoxSpacingSpinner.LowerLimitInclusive = false;
      app.VoxSpacingSpinner.Step = 0.1;
      app.VoxSpacingSpinner.ValueChangedFcn = createCallbackFcn(app, ...
        @ValueChangedCb, true);
      
      % Vox: Scale
      app.VoxScaleSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
        'label', 'Voxel Scale');
      app.VoxScaleSpinner.Layout.Row = 5;
      app.VoxScaleSpinner.Layout.Column = 1;
      app.VoxScaleSpinner.Limits = [0, Inf];
      app.VoxScaleSpinner.Step = 1;
      app.VoxScaleSpinner.LowerLimitInclusive = false;
      app.VoxScaleSpinner.ValueChangedFcn = createCallbackFcn(app, ...
        @ValueChangedCb, true);
      
      % Vox: Even range
      app.VoxEvenRange = uicheckbox(app.MainGrid);
      app.VoxEvenRange.Text = 'Voxel Even Range';
      app.VoxEvenRange.Layout.Row = 6;
      app.VoxEvenRange.Layout.Column = 1;
      app.VoxEvenRange.ValueChangedFcn = createCallbackFcn(app, ...
        @ValueChangedCb, true);
      
      % Surf: show normals
      app.SurfNormals = uicheckbox(app.MainGrid);
      app.SurfNormals.Text = 'Show Normals';
      app.SurfNormals.Layout.Row = 7;
      app.SurfNormals.Layout.Column = 1;
      app.SurfNormals.ValueChangedFcn = createCallbackFcn(app, ...
        @SurfNormalsValueChangedCb, true);
      
      % Surf: Normals scale
      app.SurfNormalScale = ott.ui.support.LabeledSpinner(app.MainGrid, ...
        'label', 'Normal scale');
      app.SurfNormalScale.Layout.Row = 8;
      app.SurfNormalScale.Layout.Column = 1;
      app.SurfNormalScale.Limits = [0, Inf];
      app.SurfNormalScale.LowerLimitInclusive = false;
      app.SurfNormalScale.Step = 0.1;
      app.SurfNormalScale.ValueChangedFcn = createCallbackFcn(app, ...
        @ValueChangedCb, true);
      
      % Surf: show edges
      app.SurfEdges = uicheckbox(app.MainGrid);
      app.SurfEdges.Text = 'Show Edges';
      app.SurfEdges.Layout.Row = 9;
      app.SurfEdges.Layout.Column = 1;
      app.SurfEdges.ValueChangedFcn = createCallbackFcn(app, ...
        @ValueChangedCb, true);
      
      % Surf: shading smooth
      app.SurfSmooth = uicheckbox(app.MainGrid);
      app.SurfSmooth.Text = 'Smooth Shading';
      app.SurfSmooth.Layout.Row = 10;
      app.SurfSmooth.Layout.Column = 1;
      app.SurfSmooth.ValueChangedFcn = createCallbackFcn(app, ...
        @ValueChangedCb, true);
      
      % Update button
      app.UpdateButton = ott.ui.support.UpdateCheckButton(app.MainGrid);
      app.UpdateButton.Layout.Row = 12;
      app.UpdateButton.Layout.Column = 1;
      addlistener(app.UpdateButton, "UpdateCalled", ...
          @(h,e) app.UpdateCb(e));
      
    end
    
    function createRightComponents(app)

      % Create UIAxes
      app.UIAxes = uiaxes(app.RightPanel);
      xlabel(app.UIAxes, '')
      ylabel(app.UIAxes, '')
      app.UIAxes.Position = [10 10 373 app.windowSize(2)-20];
      
      % Create loading text
      app.LoadingText = uilabel(app.RightPanel);
      app.LoadingText.Position = [150, app.windowSize(2)/2, 70, 22];
      app.LoadingText.Text = 'Loading...';
      app.LoadingText.Visible = 'off';
      app.LoadingText.BackgroundColor = app.UIFigure.Color;
      app.LoadingText.HorizontalAlignment = 'center';
    end
  end

  % App creation and deletion
  methods (Access = public)

    function app = Visualise()
      % Start the shape visualisation interface
      
      app = app@ott.ui.support.AppTwoColumn();
      
      if nargout == 0
        clear app;
      end
    end
  end
end