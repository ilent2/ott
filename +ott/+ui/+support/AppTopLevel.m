classdef AppTopLevel < ott.ui.support.AppBase
% Defines the properties expected by apps used in the Launcher.
%
% Abstract properties:
%   - cnameText
%   - nameText
%   - aboutText
%   - helpText
%
% Adds a menu to the window.  If the class is also a GenerateCodeMenu
% instance, adds a GerenateCode button to the menu.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Abstract, Constant)
    cnameText
    nameText
    aboutText
    helpText
  end
  
  methods (Abstract, Access=protected)
    createMainComponents(app)
    setDefaultValues(app, evt)
  end
  
  methods (Access=protected)
    function createFileMenuContent(app)
      % Overload this method to add other file menu content
    end
  end
  
  methods (Access=protected)
    function createHelpMenu(app)
      % Add a OTT help menu to the window

      % Create HelpMenu
      menu = uimenu(app.UIFigure);
      menu.Text = 'Help';

      % Create OTTonGitHubMenu
      OTTonGitHubMenu = uimenu(menu);
      OTTonGitHubMenu.Text = 'OTT on GitHub';
      OTTonGitHubMenu.MenuSelectedFcn = ...
          @(~, ~) web('https://github.com/ilent2/ott', '-browser');

      % Create OnlineManualMenu
      OnlineManualMenu = uimenu(menu);
      OnlineManualMenu.Text = 'Online Manual';
      OnlineManualMenu.MenuSelectedFcn = ...
          @(~, ~) web('https://ott.readthedocs.io/en/latest/', '-browser');

      % Create ReportaBugMenu
      ReportaBugMenu = uimenu(menu);
      ReportaBugMenu.Text = 'Report a Bug';
      ReportaBugMenu.MenuSelectedFcn = ...
          @(~, ~) web('https://github.com/ilent2/ott/issues', '-browser');

      % Create window help popup button
      menuItem = uimenu(menu);
      menuItem.Separator = 'on';
      menuItem.Text = 'Window Help';
      menuItem.MenuSelectedFcn = @(~, ~) helpdlg(app.helpText);
    end
    
    function createFileMenu(app)
      % Create the file menu

      % Create HelpMenu
      menu = uimenu(app.UIFigure);
      menu.Text = 'File';
      
      % Create ResettoDefaultsMenu
      menuItem = uimenu(menu);
      menuItem.MenuSelectedFcn = ...
        createCallbackFcn(app, @setDefaultValues, true);
      menuItem.Text = 'Reset to Defaults';
      
      % Refresh inputs
      if isa(app, 'ott.ui.support.RefreshInputsMenu')
        menuItem = uimenu(menu);
        menuItem.MenuSelectedFcn = ...
          createCallbackFcn(app, @refreshInputsCb, true);
        menuItem.Text = 'Refresh Inputs';
      end

      % Create GenerateCodeMenu
      if isa(app, 'ott.ui.support.GenerateCodeMenu')
        menuItem = uimenu(menu);
        menuItem.Separator = 'on';
        menuItem.MenuSelectedFcn = createCallbackFcn(app, ...
            @GenerateCodeMenuSelected, true);
        menuItem.Text = 'Generate Code...';
      end
      
      % Add other content
      app.createFileMenuContent();
      
      % Add show launcher button
      menuItem = uimenu(menu);
      menuItem.Separator = 'on';
      menuItem.Text = 'Show OTT Launcher';
      menuItem.MenuSelectedFcn = @(~, ~) ott.ui.Launcher.LaunchOrPresent;
      
    end
    
    function createMenuComponents(app)
      
      % Add a menu bar to the window
      app.createFileMenu();
      app.createHelpMenu();
      
    end
    
    function createComponents(app)
      
      app.createMenuComponents();
      app.createMainComponents();
      
      % Populate drop down variable selectors
      if isa(app, 'ott.ui.support.RefreshInputsMenu')
        app.refreshInputsCb();
      end
      
      % Set default values of main components
      app.setDefaultValues();
      
    end
  end
end
