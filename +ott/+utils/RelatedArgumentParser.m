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
%   p.addRule('index_medium = speed0 ./ speed');
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
%   - Rules      -- Relationships between rules (symbolic expressions)
%   - Defaults   -- Default values for after rule evaluation
%   - RequiredResults -- Results of required arguments

  properties
    Required    % Names of required arguments (char arrays)
    Optional    % Names of optional arguments (char arrays)
    Rules       % Relationships between rules (symbolic expressions)
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
      p.Rules = {};
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

    function addRule(p, expr)
      % Add a rule for calculating a required argument
      %
      % Usage
      %   p.addRule(sym)
      %
      %   p.addRule(string)  As above, but uses ``str2sym`` to
      %   convert the string to a symbolic expression.

      % Convert from string to symbolic expression
      if isstr(expr)
        expr = str2sym(expr);
      end

      % Store the rule
      p.Rules{end+1} = expr;
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

      % Get list of variables we want to solve for
      solve_vars = {};
      for ii = 1:length(p.Required)
        rname = p.Required{ii};
        if isempty(p.RequiredResults.(rname))
          solve_vars{end+1} = rname;
        end
      end

      % Check if we have work to do
      if isempty(solve_vars)
        return;
      end

      % Get a list of equations to solve
      eqns = p.Rules;

      % Get a list of known variables (for substitution)
      % Reshape known_vals with an extra dimension for later substitution
      known_syms = {};
      known_vals = {};
      unknown_syms = {};
      for ii = 1:length(p.Required)
        rname = p.Required{ii};
        val = p.RequiredResults.(rname);
        if ~isempty(p.RequiredResults.(rname))
          known_syms{end+1} = sym(rname);
          known_vals{end+1} = reshape(val, [1, size(val)]);
        end
      end
      for ii = 1:length(p.Optional)
        rname = p.Optional{ii};
        val = p.Results.(rname);
        if ~isempty(p.Results.(rname))
          known_syms{end+1} = sym(rname);
          known_vals{end+1} = reshape(val, [1, size(val)]);
        else
          unknown_syms{end+1} = sym(rname);
        end
      end

      % Solve the equations
      WS = warning('off', 'symbolic:solve:warnmsg3');
      results = solve(eqns, [unknown_syms, solve_vars], ...
          'ReturnConditions', true);
      warning(WS);

      % Check for used parameters and assign results
      for ii = 1:length(solve_vars)

        % Get solution for this variable
        var_result = results.(solve_vars{ii});
        if isempty(var_result)
          continue;
        end

        % Apply substitutions
        var_result_subs = subs(var_result, known_syms, known_vals);
        
        % Discard non-unique solutions from solve
        sz = size(var_result_subs);
        var_result_subs = reshape(var_result_subs, sz(1), []);
        var_result_subs = unique(var_result_subs, 'rows');
        
        % Check we have a unique solution
        if size(var_result_subs, 1) > 1
          warning('ott:utils:RelatedArgumentParser:non_unique_soln', ...
            'Unable to find unique solution for parameter');
          var_result_subs = var_result_subs(1, :);
        end
        
        % Reshape to original variable size (now we have applied unique)
        % Drop first dimension which we added to known_vals
        var_result_subs = reshape(var_result_subs, [sz(2:end), 1]);
        
        % Evaluate result
        result = p.castToDouble(var_result_subs);

        % Check result
        if ~isempty(result)
          used_knowns = {};
          for jj = 1:length(known_syms)
            if any(symvar(var_result) == known_syms{jj})
              used_knowns{end+1} = char(known_syms{jj});
            end
          end
          used_parameters = [used_parameters, used_knowns];

          p.RequiredResults.(solve_vars{ii}) = result;
        end
      end
    end

    function ret = castToDouble(p, val)
      % Trys to cast the variable to a double, if it can't, returns []
      %
      % There might be a more optimal way to do this cast without
      % catching an exception, but I can't find it in the Matlab doc.

      % Try and convert the number to a double
      ret = [];
      try
        ret = double(val);
      catch ME
        switch ME.identifier
          case 'symbolic:double:cantconvert'
            % Nothing to do
          otherwise
            rethrow(ME);
        end
      end
    end
  end
end
