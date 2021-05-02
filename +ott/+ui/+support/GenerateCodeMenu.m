classdef GenerateCodeMenu < handle
% Inherit from this class if your application has a GenerateCode menu.
%
% Abstract methods
%   - generateCode

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Abstract, Access=protected)
    generateCode(app)
  end
  
  methods (Access={?ott.ui.support.AppTopLevel})
    function GenerateCodeMenuSelected(app, ~)
      
      % Get code
      mcode = app.generateCode();
      
      % Add generation info
      mcode = [{'% Generated using the optical tweezers toolbox.'}, {''}, mcode];
      
      % Convert to char array, add line endings
      mcode_str = [];
      for ii = 1:length(mcode)
          mcode_str = [mcode_str, mcode{ii}, newline]; %#ok<AGROW>
      end

      % Create new document
      editorDoc = matlab.desktop.editor.newDocument(mcode_str); 
      editorDoc.smartIndentContents();
      editorDoc.goToLine(1);
      
    end
  end
end