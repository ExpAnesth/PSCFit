% Template script calling function sumPlot_PSC which generates simple
% collections of line plots of PSC parameters, one figure per dataSet

% ---- paths, printing prefs
compName=lower(getenv('computername'));
switch compName
  case {'hh-i7'}
    dataPath='e:\_data\_MatlabMakeover\ipsc\results_figures\';
    plotPath='e:\_data\_MatlabMakeover\ipsc\results_figures\';
  case {'hh64'}
    dataPath='e:\_data\_MatlabMakeover\ipsc\results_figures\';
    plotPath='e:\_data\_MatlabMakeover\ipsc\results_figures\';
  otherwise
    error('unknown machine')
end

% data set; column order:
% 1. data file name
% 2. indexes into r.pscrMn, corresponding to drug conditions listed in
%    ds.indepPar and ds.indepParLabel, to be compared
% 3. axis labels
dataSet={...
  '\substanceX.mat', [1 2], {'ctrl', 'drug X (0.25 µM)'};...
  '\substanceY.mat', [1 2], {'ctrl', 'drug Y (0.5 µM)'};...  
  };

% parameters to plot; column order: parameter name, axis label, log or lin
% scale (see template_call_pscdeal for a full list of parameters that are
% available)
fullPSCPar={...
  'allAmp', 'peak amplitude (pA)',         'lin';
  'tDecay', '{\tau}_{decay} (ms)',         'lin';
  'chargePPsc', 'charge per PSC (fC)',  'lin';
  'chargePhas', 'charge/time (pA)', 'log';
  };
nPar=size(fullPSCPar,1);

% some parameters governing plot appearance and export
ap.subPlotGrid=[3 2];
% print?
ap.printas='-dpsc2';
ap.printas=[];

% graphics scaling
labelscale('fontSz',10,'scaleFac',1,'lineW',.5,'markSz',7);

% call 
sumPlot_PSC(dataPath,plotPath,dataSet,fullPSCPar,ap);