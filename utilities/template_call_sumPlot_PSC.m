% Template script calling function sumPlot_PSC which generates simple
% collections of line plots of PSC parameters, one figure per dataSet

% ---- paths, printing prefs
compName=lower(getenv('computername'));
switch compName
  case {'hh-i7'}
    dataPath='d:\hh\projects\ctx_propoSevo\data_IPSC\';
    plotPath='d:\hh\projects\ctx_propoSevo\rawFig\';
  case {'hh64'}
    dataPath='d:\hh\projects\ctx_propoSevo\data_IPSC\';
    plotPath='d:\hh\projects\ctx_propoSevo\rawFig\';
  case {'eval-lmb'}
    dataPath='d:\hh\projects\ctx_propoSevo\data_IPSC\';
    plotPath='d:\hh\projects\ctx_propoSevo\rawFig\';
  otherwise
    error('unknown machine')
end


% data set; column order:
% 1. data file name
% 2. indexes into PSCRMN, corresponding to drug conditions listed in
%    ds.indepPar and ds.indepParLabel, to be compared (second will be
%    normalized by first)
% 3. axis labels
dataSet={...
  '\IPSC_prop.mat', [2 4 5], {'ctrl','0.5 µm prop.','wash'};...
  '\IPSC_sevo.mat', [2 4 5], {'ctrl','70 µM sevo','wash'};...
  '\IPSC_sevoprop.mat', [2 4 5], {'ctrl','70 µM sevo+0.5 µm prop.','wash'};...
  };


% parameters to plot, axis labels, log or lin, detected or fitted

% 'tRise'    'tDecay'    'width'    'amp'    'allTsl'    'tsl'    'freq'
% 'freqFit'    'chargePPsc' 'chargePscTot'    'baseline'    'noise'
% 'chargePhas'    'allAmp'    'allTRise20_80'    'tRise20_80'

fullPSCPar={...
  'freq',          'frequency (Hz)',       'log', 'detected';
  %   'freqFit',       'frequency (Hz)',       'log', 'fitted';
  'allAmp', 'peak amplitude (pA)',         'lin', 'detected';
  %   'amp', 'peak amplitude (pA)',            'lin', 'fitted';
  % 'allTRise20_80', '20-80% rise time (ms)','lin', 'detected';
  %   'tRise20_80' , '20-80% rise time (ms)',  'lin', 'fitted' ;
  'tDecay', '{\tau}_{decay} (ms)',         'lin', 'fitted';
  'chargePPsc', 'charge per PSC (fC)',  'lin', 'fitted';
  'chargePhas', 'charge/time (pA)', 'log', 'detected';
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