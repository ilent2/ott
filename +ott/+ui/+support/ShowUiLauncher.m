function ShowUiLauncher()
% Looks for an active UI launcher and shows it, otherwise opens a new one

  % Search for a figure with the specific name
  LauncherName = 'Optical Tweezers Toolbox Launcher';
  lh = findall(0, 'Name', LauncherName);
  
  if isempty(lh)
    % Open a new UI Launcher
    ott.ui.Launcher();
  else
    % Bring the launcher to the forground
    figure(lh(1));
  end

end