function menu = createHelpMenu(fig)
% Add a OTT help menu to the specified uifigure
%
% Usage
%   menu = createHelpMenu(fig)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% Create HelpMenu
menu = uimenu(fig);
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

end