% Run OTT tests

addpath('../../ott');

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin

showCoverage = false;

% Array of sub-packages to test
subpackages = {};
subpackages{end+1} = 'beam';
subpackages{end+1} = 'bsc';
subpackages{end+1} = 'drag';
subpackages{end+1} = 'particle';
subpackages{end+1} = 'shape';
subpackages{end+1} = 'tmatrix';
subpackages{end+1} = 'tmatrix/dda';
subpackages{end+1} = 'tools';
% subpackages{end+1} = 'ui';
subpackages{end+1} = 'utils';
% subpackages{end+1} = 'examples';

% Array of sub-packages to generate coverage reports for
coverage = {};
if showCoverage
  coverage{end+1} = 'ott.beam';
  coverage{end+1} = 'ott.bsc';
  coverage{end+1} = 'ott.drag';
  coverage{end+1} = 'ott.particle';
  coverage{end+1} = 'ott.shape';
  coverage{end+1} = 'ott.tmatrix';
  coverage{end+1} = 'ott.tools';
  %coverage{end+1} = 'ott.ui';
  coverage{end+1} = 'ott.utils';
end

% Generate coverage report runner for packages
runner = TestRunner.withTextOutput;
if ~isempty(coverage)
  runner.addPlugin(CodeCoveragePlugin.forPackage(...
      coverage, 'IncludingSubpackages', true));
end

if showCoverage
  % Add coverage report for examples directory
  runner.addPlugin(CodeCoveragePlugin.forFolder('../examples', ...
      'IncludingSubfolders', true));
end

% Build test suite
suites = {};
for ii = 1:numel(subpackages)
  suites{ii} = TestSuite.fromFolder(pwd, 'IncludingSubfolders', true, ...
      'BaseFolder', fullfile(pwd, subpackages{ii}));
end

% Generate coverage report
result = runner.run(horzcat(suites{:}));

%% Display some stats

rt = table(result);
sortrows(rt,'Duration')

disp(result)
