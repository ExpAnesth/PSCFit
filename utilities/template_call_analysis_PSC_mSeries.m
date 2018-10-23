% Template script calling function analysis_PSC_mSeries which generates a
% boxplot of normalized PSC parameter values for one or several data
% series; the script also computes statistics on the differences between
% the series

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
%    ds.indepPar and ds.indepParLabel, to be compared (second will be
%    normalized by first)
% 3. label (to be used for legend)
% 4. plot color 
% 5. anything else (currently, the shortcut needed for looking up
%    properties of the data set like associated color)
dataSet={...
  '\substanceX.mat', [1 2], 'drug X (0.25 µM)', [0 0 1], '';...
  '\substanceY.mat', [1 2], 'drug Y (0.5 µM)',  [1 0 1], '';...  
  };

% parameters to plot; column order: parameter name, axis label (see
% template_call_pscdeal for a full list of parameters that are available)
fullPSCPar={...
  'allAmp', 'peak ampl.';
  'tDecay', 'decay time';
  'chargePPsc', 'charge/PSC';
  };
nPar=size(fullPSCPar,1);

% some parameters governing box plot appearance and export
ap.subPlotGrid=[2 2];
ap.dataLim=[0 3.3];
ap.factorGap=[7 2];
ap.whisker=1.5;
ap.plotFn='boxplot_IPSCPar';
% print?
ap.printas='-dpsc2';
ap.printas=[];

% graphics scaling
labelscale('fontSz',9,'scaleFac',1,'lineW',.25,'markSz',3.5);

%% call to analysis_PSC_mSeries
mkfig(1,'b');
clf
[dCell1way,gCell1way,dCell2way,gCell2way]=...
  analysis_PSC_mSeries(dataPath,plotPath,dataSet,fullPSCPar,ap);


%% oneway stats 
for pIx=1:nPar
  disp([fullPSCPar{pIx,1} ', one-way analysis: '])
  stats=mes1way(dCell1way{pIx},{'partialeta2','g_psi'},'group',gCell1way{pIx},'cWeight',[-1 1 0]);
end

%% twoway stats
for pIx=4% 1:nPar
  disp([fullPSCPar{pIx,1} ', two-way analysis: '])
  cw=[1 -1 0; -1 1 0]'; % sevo vs prop
%   cw=[1 0 -1; -1 0 1]'; % sevo vs prop+sevo  
%   cw=[0 1 -1; 0 -1 1]'; % prop vs prop+sevo  
  stats=mes2way(dCell2way{pIx},gCell2way{pIx},{'partialeta2','g_psi'},...
    'fName',{'series','drug'},'isDep',[0 1],'cWeight',cw,'nBoot',5000,'doDataPlot',true);
  
end
  
  
  
