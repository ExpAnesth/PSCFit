%% -------------------- PREPARATIONS --------------------------------------
compName=lower(getenv('computername'));
switch compName
  case {'hh-i7'}
    ds.dataPath='d:\_data\otc_ctx\ACh\AChBlockIPSC\';
    ds.plotPath='d:\hh\projects\ctx_ACh\rawFig\';
  case {'hh64','hh-i5'}
    ds.dataPath='e:\_data\otc_ctx\ACh\AChBlockIPSC\';
    ds.plotPath='d:\hh\projects\ctx_ACh\rawFig\';
  case {'eval-lmb'}
    ds.dataPath='h:\_data\otc_ctx\ACh\AChBlockIPSC\';
    ds.plotPath='d:\hh\projects\ctx_ACh\rawFig\';
  otherwise
    error('machine not defined');
end

% --- graphics settings
labelscale('fontSz',7,'scaleFac',1,'lineW',.3,'markSz',3);
ds.printas='-djpeg98';
ds.printas='-dpsc2';
% ds.printas=[];
% color and symbols for different ACh conditions
fac=graphicsDef_ACh;
% scaling factor for bubble plots: an amplitude of 1000 pA shall correspond
% to area 100;
ds.ampScalFac=1/10;
% if true, surface/contour plots will have same color scaling across drug
% conditions
ds.doSameColorScale=false;
% if true, all surface/contour/scatter plots will have axis limits
% embracing the data (as opposed to limits imposed by ds.fullPSCPar below)
ds.doIndividualAxLim=true;

% --- analysis parameters
% - parameters to load, check and possibly analyze
% - description including units
% - bins for 2D hist (at the same time, these define the limits of
%  acceptable values)
% - ticks for axes
% - type of par: - 'detected': from all detected IPSCs, including those which could be fitted
%                - 'fitted': from fitted IPSCs only
%                          *** DETECTED MUST BE LISTED FIRST ***
ds.fullPSCPar={...
  'allTRise20_80' , '20-80% rise time (ms)',2.^[log2(.05):.22142:log2(5)],[.1 1], 'detected';
  'allAmp', 'peak amplitude (pA)',2.^[log2(10):.2214614:log2(1000)],[10 100 1000],'detected';
  'tDecay', '{\tau}_{decay} (ms)',2.^[log2(2):.188128:log2(100)],[10 100],        'fitted';
  'amp', 'peak amplitude (pA)',2.^[log2(10):.2214614:log2(1000)],[10 100 1000],   'fitted';
  'tRise20_80' , '20-80% rise time (ms)',2.^[log2(.05):.22142:log2(5)],[.1 1],    'fitted' ;
  'chargePPsc', 'charge per PSC (pA*ms)',2.^[log2(0.1):.21:log2(8)],[.1 1],       'fitted'; ... % §§ check limits!
};

% set to zero for no bootstrapping/randomization 
ds.nBoot=1000;

% ** histogram analysis parameters
% - percentiles of null distribution of bin-wise differences to mark in
% plots (values, color)
ds.ndPrctile={2.5,[0 .8 0];97.5,[.9 .9 .2]};


%%
% subdir
ds.dSubDir='Ctrl-ACh-Block\Figs\';
% name of data file
ds.dFn='AChBlock.mat';
% indexes into PSCRMN, corresponding to drug conditions listed in
% ds.indepPar and ds.indepParLabel, to be analyzed
ds.drugIx=[1 2 3];
ds.drugTag={'ACh0','ACh+','ACh-'};
% indexes into PSCRMN to drug conditions to be used for statistical
% comparisons (must be a subset of the above; order matters)
ds.drugStatsIx=[2 1];
% *** index to two PARAMETERS (as listed in ds.fullPSCPar) to be used for
% primary analysis (construction of 2D histograms & identification of
% altered spots in parameter space or clustering)
ds.anPSCParInd=[1 2];
% index to parameters to use in scatter plots with difference scores on
% abscissa
ds.anPSCParInd2=[3 6];

% ** cluster analysis parameters:
% - do it?
ds.doClusterAnalysis=true;
% - indexes into PSCRMN to drug condition(s) to be used for cluster detection/generation
ds.drugClusterIx=[1 2 3];
% - shall data be log-transformed for clustering? 
ds.doLogTransform=true;

% define colors corresponding to ACh status 
fIx=find(strcmp('ACh status',{fac.name}));
ds.pCol=[];
for dcIx=1:numel(ds.drugIx)
  ds.pCol(dcIx,:)=fac(fIx).color(strcmp(ds.drugTag{dcIx},fac(fIx).levelName),:);
end

% *** call analysis function
analysis_IPSC_mvar(ds);

return

%% 
% subdir
ds.dSubDir='ACh_ACh_Block\Figs\';
% name of data file
ds.dFn='AChAChBlock.mat';
% indexes into PSCRMN, corresponding to drug conditions listed in
% ds.indepPar and ds.indepParLabel, to be analyzed
ds.drugIx=[1 2 3];
ds.drugTag={'ACh0','ACh+','ACh-'};
% indexes into PSCRMN to drug conditions to be used for statistical
% comparisons (must be a subset of the above; order matters)
ds.drugStatsIx=[2 3];
% *** index to two PARAMETERS (as listed in ds.fullPSCPar) to be used for
% primary analysis (construction of 2D histograms & identification of
% altered spots in parameter space or clustering)
ds.anPSCParInd=[1 2];
% index to parameters to use in scatter plots with difference scores on
% abscissa
ds.anPSCParInd2=[1 2 3 6];

% cluster analysis parameters:
% - do it?
ds.doClusterAnalysis=false;
% - indexes into PSCRMN to drug condition(s) to be used for cluster detection/generation
ds.drugClusterIx=[1 2 3];
% - shall data be log-transformed for clustering? 
ds.doLogTransform=true;

% define colors corresponding to ACh status 
fIx=find(strcmp('ACh status',{fac.name}));
ds.pCol=[];
for dcIx=1:numel(ds.drugIx)
  ds.pCol(dcIx,:)=fac(fIx).color(strcmp(ds.drugTag{dcIx},fac(fIx).levelName),:);
end

% *** call analysis function
analysis_IPSC_mvar(ds);

%%




% % subdir
% ds.dSubDir='ACh-Dia-Block\Figs\';
% % name of data file
% ds.dFn='AChDiaBlock.mat';
% % indexes into PSCRMN, corresponding to drug conditions listed in
% % ds.indepPar and ds.indepParLabel
% ds.drugIx=[1 2 3];
% ds.drugTag={'ACh0','ACh+','ACh-'};
%  
% % subdir
% ds.dSubDir='ACh-Zolpi200-Block\Figs\';
% % name of data file
% ds.dFn='AChZolpi200Block.mat';
% % indexes into PSCRMN, corresponding to drug conditions listed in
% % ds.indepPar and ds.indepParLabel
% ds.drugIx=[1 2 3];
% ds.drugTag={'ACh0','ACh+','ACh-'};
% 
% % subdir
% ds.dSubDir='ACh_TTX\Figs\';
% % name of data file
% ds.dFn='AChTTX.mat';
% % indexes into PSCRMN, corresponding to drug conditions listed in
% % ds.indepPar and ds.indepParLabel
% ds.drugIx=[1 2 3];
% ds.drugTag={'ACh0','ACh+','ACh-'};


% % subdir
% ds.dSubDir='ACh-Zolpi-Block\Figs\';
% % name of data file
% ds.dFn='AChZolpiBlock.mat';
% % indexes into PSCRMN, corresponding to drug conditions listed in
% % ds.indepPar and ds.indepParLabel
% ds.drugIx=[1 2 3];

