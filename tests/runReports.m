% This file runs the code analyser reports for each toolbox package

reportlist = {};
reportlist{end+1} = '../+ott';
reportlist{end+1} = '../+ott/+beam';
reportlist{end+1} = '../+ott/+bsc';
reportlist{end+1} = '../+ott/+particle';
reportlist{end+1} = '../+ott/+tmatrix';
reportlist{end+1} = '../+ott/+shape';
reportlist{end+1} = '../+ott/+shape/+mixin';
reportlist{end+1} = '../+ott/+drag';
reportlist{end+1} = '../+ott/+tools';
reportlist{end+1} = '../+ott/+ui';
reportlist{end+1} = '../+ott/+utils';
reportlist{end+1} = '../examples';

for ii = 1:numel(reportlist)
  web('', '-new');
  mlintrpt(reportlist{ii},'dir');

  web('', '-new');
  contentsrpt(reportlist{ii});
end
