classdef GridWidget < handle
% A container like class with a top-level grid

  properties
    Grid          matlab.ui.container.GridLayout
  end
  
  properties (Dependent)
    Layout
  end
  
  methods
    function obj = GridWidget(parent, varargin)
      
      p = inputParser;
      p.parse(varargin{:});
      
      % Create grid
      obj.Grid = uigridlayout(parent);
      obj.Grid.Padding = [0 0 0 0];
      
    end
  end
  
  methods
    function val = get.Layout(obj)
      val = obj.Grid.Layout;
    end
    
    function set.Layout(obj, val)
      obj.Grid.Layout = val;
    end
  end
end