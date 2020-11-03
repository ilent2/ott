function menuItem = createLauncherMenuItem(parent)
% Add a Show Launcher menu item to the specified menu
%
% Usage
%   menuItem = createLauncherMenuItem(parent)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% Create menu item
menuItem = uimenu(parent);
menuItem.Separator = 'on';
menuItem.Text = 'Show OTT Launcher';
menuItem.MenuSelectedFcn = @ott.ui.Launcher.LaunchOrPresent;

end

