% Run OTT tests

addpath('../../ott');

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin

% Array of sub-packages to test
subpackages = {};
% subpackages{end+1} = 'bsc';
subpackages{end+1} = 'tmatrix';
subpackages{end+1} = 'tmatrix/dda';
% subpackages{end+1} = 'shape';
% subpackages{end+1} = 'drag';
% subpackages{end+1} = 'utils';
% subpackages{end+1} = 'ui';
% subpackages{end+1} = 'examples';

% Array of sub-packages to generate coverage reports for
coverage = {};
% coverage{end+1} = 'ott.bsc';
coverage{end+1} = 'ott.tmatrix';
% coverage{end+1} = 'ott.shape';
% coverage{end+1} = 'ott.drag';
% coverage{end+1} = 'ott.utils';
%coverage{end+1} = 'ott.ui';
%coverage{end+1} = 'examples';

% Generate coverage report runner
runner = TestRunner.withTextOutput;
if ~isempty(coverage)
  runner.addPlugin(CodeCoveragePlugin.forPackage(...
      coverage, 'IncludingSubpackages', true))
end

% Build test suite
suites = {};
for ii = 1:numel(subpackages)
  suites{ii} = TestSuite.fromFolder(pwd, 'IncludingSubfolders', true, ...
      'BaseFolder', fullfile(pwd, subpackages{ii}));
end

% Generate coverage report
result = runner.run(horzcat(suites{:}))

