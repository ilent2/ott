classdef RelatedArgumentParser < inputParser
% A class for parsing arguments related by rules
%
% This class is a type of input parser which attempts to provide
% the required arguments by apply a set of rules to the required
% and optional arguments.  This can be useful when constructing,
% for instance, a beam where the refractive index can be specified
% either directly, as a ratio of medium speed and vacuum speed,
% or as a permittivity.
%
% Example
% .. code:: matlab
%
%   p = ott.utils.RelatedArgumentParser;
%   p.addRequired('index_medium', 1.0);
%   p.addOptional('speed0', 2.0);
%   p.addOptional('speed');
%   p.addRule('index_medium', @(s1, s0) s0 ./ s1, 'speed', 'speed0');
%   p.parse('speed', 1.0);
%   disp(p.Results.index_medium)  % should display 2
%
% Methods
%   - addRequired  -- Required named argument
%   - addOptional  -- Optional named argument
%   - addRule      -- Add a rule relating arguments
%   - parse        -- Parse the arguments and apply rules
%
% Properties
%   - Required   -- Names of required arguments (char arrays)
%   - Optional   -- Names of optional arguments (char arrays)
%   - Rules      -- Relationships between rules (struct of function_handle's)
%   - Defaults   -- Default values for after rule evaluation
%   - RequiredResults -- Results of required arguments

  properties
    Required    % Names of required arguments (char arrays)
    Optional    % Names of optional arguments (char arrays)
    Rules       % Relationships between rules (function_handle's)
    Defaults    % Default values for after rule evaluation
    RequiredResults % Results of required arguments
  end

  methods
    function p = RelatedArgumentParser()
      % Construct a new related argument parser

      % Call super-class
      p = p@inputParser();

      % Initialize properties
      p.Required = {};
      p.Optional = {};
      p.Rules = struct();
      p.Defaults = struct();
      p.RequiredResults = struct();
    end

    function addRequired(p, name, default_value, varargin)
      % Add a required parameter to the input parser
      %
      % Usage
      %   p.addRequired(name)
      %
      %   p.addRequired(name, default_value, ...)
      %   Specifies a default value to assign after applying relations.
      %   Additional arguments are passed to addParameter

      % Call base method to add parameter
      p.addParameter(name, [], varargin{:});

      % Store name
      p.Required{end+1} = name;

      % Setup rules cell array
      p.Rules.(name) = {};

      % Store post-relationship defaults
      if nargin >= 3
        p.Defaults.(name) = default_value;
      end
    end

    function addOptional(p, name, default_value, varargin)
      % Add an optional named argument to the input parser
      %
      % Usage
      %   p.addOptional(name)
      %
      %   p.addOptional(name, default_value, ...)
      %   Specifies the default value to assign before applying relations.
      %   Additional arguments are passed to addParameter
      
      % Handle default value
      if nargin < 3
        default_value = [];
      end
      
      % Store name
      p.Optional{end+1} = name;

      % Call base method to add parameter
      p.addParameter(name, default_value, varargin{:});
    end

    function addRule(p, var_name, func, varargin)
      % Add a rule for calculating a required argument
      %
      % Usage
      %   p.addRule(var_name, fun, [arg1, arg2, ...])
      %
      % Parameters
      %   - var_name -- Name of the required argument
      %   - func -- Function handle describing rule
      %   - arg1, arg2, ... -- Named arguments to pass to function

      % Verify var_name is a required argument
      assert(any(strcmpi(var_name, p.Required)), ...
        'ott:utils:RelatedArgumentParser:var_not_a_required', ...
        'var_name must be a required argument');

      % Verify arguments are parameters
      for ii = 1:length(varargin)
        assert(any(strcmpi(varargin{ii}, p.Parameters)), ...
          'ott:utils:RelatedArgumentParser:arg_not_a_parameter', ...
          ['arg ' varargin{ii} ' not a valid parameter name']);
      end

      % Store rule
      p.Rules.(var_name){end+1} = {func, varargin};
    end

    function parse(p, varargin)
      % Parse arguments and apply relationship rules

      % Start by parsing the arguments
      parse@inputParser(p, varargin{:});

      % Set-up required results array
      for ii = 1:length(p.Required)
        p.RequiredResults.(p.Required{ii}) = p.Results.(p.Required{ii});
      end

      % Apply the relationship rules
      used_parameters = p.resolve({});

      % Set empty required arguments to default value
      fnames = fieldnames(p.Defaults);
      for ii = 1:length(fnames)
        if isempty(p.RequiredResults.(fnames{ii}))
          p.RequiredResults.(fnames{ii}) = p.Defaults.(fnames{ii});
        end
      end

      % Apply the relationship rules (to catch any dependent on defaults)
      used_parameters = p.resolve(used_parameters);

      % Check for extra arguments and warn
      for ii = 1:length(p.Optional)

        pname = p.Optional{ii};

        % Don't care about optional arguments with default values
        if any(strcmpi(pname, p.UsingDefaults))
          continue;
        end

        % Don't care about optional arguments that were used
        if any(strcmpi(pname, used_parameters))
          continue;
        end

        warning('ott:utils:RelatedArgumentParser:not_used', ...
          ['Parameter ' pname ' set but not used in parsing']);
      end
    end
  end

  methods (Access=protected)
    function used_parameters = resolve(p, used_parameters)
      % Apply the relationship rules
      %
      % Usage
      %   used_params = p.resolve(used_params)

      while true

        % Keep track of if we did work in this loop
        done_work = false;

        % Loop over unresolved required arguments
        for ii = 1:length(p.Required)

          rname = p.Required{ii};

          % Check if we are resolved
          if ~isempty(p.RequiredResults.(rname))
            continue;
          end

          % Apply rules
          for rule = p.Rules.(rname)

            % Get function name and argument names
            func = rule{1}{1};
            arg_names = rule{1}{2};

            % Get list of arguments (from Results or RequiredResults)
            args = cell(size(arg_names));
            for jj = 1:length(arg_names)
              if any(strcmpi(arg_names{jj}, p.Required))
                args{jj} = p.RequiredResults.(arg_names{jj});
              else
                args{jj} = p.Results.(arg_names{jj});
              end
            end

            % Check we can do work
            can_work = all(cellfun(@(arg) ~isempty(arg), args));
            if ~can_work
              continue;
            end

            % Evaluate rule
            p.RequiredResults.(rname) = func(args{:});

            % Update list of optional arguments used
            used_parameters = [used_parameters, arg_names]; %#ok<AGROW>

            % Update done_work if work was done
            done_work = true;

          end
        end

        % Check if we did any work in the last iteration
        if ~done_work
          break;
        end
      end
    end
  end
end
