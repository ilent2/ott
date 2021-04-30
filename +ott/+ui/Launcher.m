classdef Launcher < matlab.apps.AppBase
% Displays a list of graphical user interfaces included with OTT.
%
% Usage
%   ott.ui.Launcher() -- Creates a new instance.
%
%   ott.ui.Launcher.LaunchOrPresent() -- Creates a new instance of presents
%   an existing instance.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    UiName = 'Optical Tweezers Toolbox Launcher';
  end

  properties (Access=public)
    LauncherUiFigure    matlab.ui.Figure
    OttOverview         matlab.ui.control.TextArea
    GitHubButton        matlab.ui.control.Button
    ManualButton        matlab.ui.control.Button
    BugButton           matlab.ui.control.Button
    LaunchButton        matlab.ui.control.Button
    TitleLabel          matlab.ui.control.Label
    ItemTextArea        matlab.ui.control.TextArea
    CategoryList        matlab.ui.control.ListBox
    ApplicationList     matlab.ui.control.ListBox
  end

  methods (Static)
    function LaunchOrPresent()
      % Launch a new Launcher instance or present an existing one

      % Look for existing instance
      lh = findall(0, 'Name', ott.ui.Launcher.UiName);

      if isempty(lh)
        ott.ui.Launcher();
      else
        figure(lh(1));
      end
    end
  end

  methods (Access=private)
    
    function ref = getCurrentAppList(app)
      
      % Meta-classes seem to need to be alocated dynamically.
      
      switch app.CategoryList.Value
        case 'Beam'
          ref = [?ott.ui.beam.PmParaxial, ...
            ?ott.ui.beam.Scattered, ?ott.ui.beam.Simple, ...
            ?ott.ui.beam.Visualise];
        case 'Drag'
          ref = [?ott.ui.drag.Simple];
        case 'Dynamics'
          ref = [?ott.ui.dynamics.Isolated];
        case 'Particle'
          ref = [?ott.ui.particle.Simple];
        case 'Shape'
          ref = [?ott.ui.shape.CadFileLoader, ?ott.ui.shape.Simple, ...
              ?ott.ui.shape.Visualise];
        case 'Tmatrix'
          ref = [?ott.ui.tmatrix.Dda, ?ott.ui.tmatrix.Mie, ...
              ?ott.ui.tmatrix.Pointmatch, ?ott.ui.tmatrix.Simple];
        case 'Tools'
          ref = [?ott.ui.tools.ForcePosition, ?ott.ui.tools.PowerSpectrum];
        otherwise
          error('Internal error');
      end
    end
    
    function LaunchButtonPushed(app, ~)
      meta = app.ApplicationList.Value;
      eval(meta.Name);
    end

    function ApplicationChangedCallback(app, ~)
      meta = app.ApplicationList.Value;
      name = findobj(meta.PropertyList, 'Name', 'nameText');
      about = findobj(meta.PropertyList, 'Name', 'aboutText');
      
      app.TitleLabel.Text = name.DefaultValue;
      app.ItemTextArea.Value = about.DefaultValue;
    end
    
    function CategoryChangedCallback(app, ~)
      ref = app.getCurrentAppList();
      
      names = {};
      for ii = 1:length(ref)
        nameProp = findobj(ref(ii).PropertyList, 'Name', 'cnameText');
        names{ii} = nameProp.DefaultValue;
      end
      
      app.ApplicationList.Items = names;
      app.ApplicationList.ItemsData = ref;
      app.ApplicationChangedCallback();
    end

    function createComponents(app)
      % Create app components (mostly based on UI designer code)

      % Create figure
      app.LauncherUiFigure = uifigure('Visible', 'off');
      app.LauncherUiFigure.Units = 'pixels';
      app.LauncherUiFigure.Position = [ott.ui.support.defaultXy, 640, 360];
      app.LauncherUiFigure.Name = ott.ui.Launcher.UiName;
      app.LauncherUiFigure.Resize = 'off';

      % Create overview box
      app.OttOverview = uitextarea(app.LauncherUiFigure);
      app.OttOverview.Editable = 'off';
      app.OttOverview.Position = [14 217 502 120];
      app.OttOverview.Value = {
        'Welcome to the Optical Tweezers Toolbox (OTT)'
        ''
        ['This is the launcher for the OTT graphical user interface.', ...
         '  This launcher provides a list of the user interfaces ', ...
         'currently included with OTT.  These interfaces are intended', ...
         ' to demonstrate the basic functionality of OTT.', ...
         '  For more details, see the online manual.']
        ''
        ['This interface and related code, except where ', ...
         'otherwise noted, is released under the Creative Commons ', ...
         'Attribution-NonCommercial 4.0 International License.']
        ''
        'If you find this package useful, please cite it as:'
        ''
        ['I. C. D. Lenton, A. B. Stilgoe, T. A. Nieminen, ', ...
         'H. Rubinsztein-Dunlop, "An object orientated optical ', ...
         'tweezers toolbox", [Journal to be decided](link to the article)']};

      % Create OTTonGitHubButton
      app.GitHubButton = uibutton(app.LauncherUiFigure, 'push');
      app.GitHubButton.Position = [528 308 100 22];
      app.GitHubButton.Text = 'OTT on GitHub';
      app.GitHubButton.ButtonPushedFcn = ...
          @(~, ~) web('https://github.com/ilent2/ott', '-browser');

      % Create OnlineManualButton
      app.ManualButton = uibutton(app.LauncherUiFigure, 'push');
      app.ManualButton.Position = [528 266 100 22];
      app.ManualButton.Text = 'Online Manual';
      app.ManualButton.ButtonPushedFcn = ...
          @(~, ~) web('https://ott.readthedocs.io/en/latest/', '-browser');

      % Create ReportaBugButton
      app.BugButton = uibutton(app.LauncherUiFigure, 'push');
      app.BugButton.Position = [528 224 100 22];
      app.BugButton.Text = 'Report a Bug';
      app.BugButton.ButtonPushedFcn = ...
          @(~, ~) web('https://github.com/ilent2/ott/issues', '-browser');

      % Create application list box
      app.ApplicationList = app.createNamedListBox(app.LauncherUiFigure, ...
          'Title', 'Application', 'Position', [140, 18, 110, 183]);
      app.ApplicationList.ValueChangedFcn = createCallbackFcn(app, ...
          @ApplicationChangedCallback, true);

      % Create category list box
      app.CategoryList = app.createNamedListBox(app.LauncherUiFigure, ...
          'Title', 'Category', 'Position', [14, 18, 110, 183]);
      app.CategoryList.Items = {'Beam', 'Drag', 'Dynamics', 'Particle', ...
         'Shape', 'Tmatrix', 'Tools'};
      app.CategoryList.ValueChangedFcn = createCallbackFcn(app, ...
          @CategoryChangedCallback, true);

      % Create Panel
      LaunchPanel = uipanel(app.LauncherUiFigure);
      LaunchPanel.Position = [267 19 361 181];
      m = 8;

      % Create LaunchButton
      app.LaunchButton = uibutton(LaunchPanel, 'push');
      app.LaunchButton.Position = [(261-m) 152 100 22];
      app.LaunchButton.Text = 'Launch';
      app.LaunchButton.ButtonPushedFcn = ...
          createCallbackFcn(app, @LaunchButtonPushed, true);

      % Create TitleLabel
      app.TitleLabel = uilabel(LaunchPanel);
      app.TitleLabel.Position = [m 152 261 22];

      % Create ItemTextArea
      app.ItemTextArea = uitextarea(LaunchPanel);
      app.ItemTextArea.Position = [m 9 (361-2*m) 137];
      app.ItemTextArea.Editable = 'off';

      % Show the figure after all components are created
      app.LauncherUiFigure.Visible = 'on';
    end

    function listBox = createNamedListBox(app, parent, varargin)
      % Create a named list box (using a uipanel and listbox

      p = inputParser;
      p.addParameter('Title', []);
      p.addParameter('Position', []);
      p.parse(varargin{:});

      % Create panel
      panel = uipanel(parent);
      panel.Title = p.Results.Title;
      panel.Position = p.Results.Position;

      % Create list box
      listBox = uilistbox(panel);
      listBox.Position = [0, 0, panel.Position(3), panel.Position(4)-20];
    end
  end

  methods (Access=public)
    function app = Launcher()
      % Create new Launcher instance or present existing instance

      % Create UI
      app.createComponents();
      
      % Run selection callbacks
      app.CategoryChangedCallback();
      
      % Register the app with App Designer
      registerApp(app, app.LauncherUiFigure)

      if nargout == 0
        clear app;
      end
    end

    % Code that executes before app deletion
    function delete(app)

      % Delete UIFigure when app is deleted
      delete(app.LauncherUiFigure);
    end
  end
end
