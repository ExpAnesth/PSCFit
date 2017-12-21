% this is a template of the input structure into pscdeal.m
% -------------------- PREPARATIONS ---------------------------------------
compName=lower(getenv('computername'));
switch compName
  case {'hh-i7'}
    dataPath='d:\_data\otc_ctx\ACh\AChBlockIPSC\Mist\ACh+Dia-Zolpi\';
    plotPath='d:\_data\otc_ctx\ACh\AChBlockIPSC\Mist\ACh+Dia-Zolpi\summary&Figures';
  case {'hh64','hh-i5'}
    dataPath='e:\_data\otc_ctx\ACh\AChBlockIPSC\Mist\ACh+Dia-Zolpi\';
    plotPath='e:\_data\otc_ctx\ACh\AChBlockIPSC\Mist\ACh+Dia-Zolpi\summary&Figures';
  case {'eval-lmb'}
    dataPath='h:\_data\otc_ctx\ACh\AChBlockIPSC\Mist\ACh+Dia-Zolpi\';
    plotPath='h:\_data\otc_ctx\ACh\AChBlockIPSC\Mist\ACh+Dia-Zolpi\summary&Figures';
  otherwise
    error('machine not defined');
end

% name of *.mat file into which to write results
resultFn='AChZolDia';
% --- data set descriptors
% the complete set of values of independent parameter (=concentrations)
ds.indepPar=[0 1 2 3];
% matching descriptors
ds.indepParLabel={'ACh','AChp_Zol','AChp_Dia','ACh_wash'};
% name suffix and extension of files produced by threshdetgui/PSCFitgui
ds.fnSuffix='_IPSC_res';
ds.fnSuffix='_cIPSC_res';

% --- analysis parameters
% -- value of independent parameter corresponding to control (against which
% others will be normalized; set to nan for no normalization)
ds.normIpVal=nan;
% -- value of independent parameter of condition against which others will
% be (statistically) COMPARED
ds.compIpVal=1;

% -- parameters of **fitted** PSCs as computed in pscfitgui (or which can
% be derived from them)
ds.pscFitPar={...
  'tRise',...           % 10-90 % rise time (ms)
  'tDecay',...          % decay time (ms)
  'width',...           % width (ms)
  'amp',...             % peak amplitude (same unit as in raw recording; as determined from fit)
  'aPeak',...           % peak amplitude (same unit as in raw recording; as determined in raw cutouts)
  'chargePPsc',...      % charge transferred per PSC (pC)
  'chargePscTot',...    % total charge transferred (same unit as in raw recording)
  'xIntvFitEnd',...     % right border of fit interval (ms)
  'qFit',...            % quality of fit (currently, R^2)
  'freqFit',...         % frequency (Hz; see parameter allFreq below)
  'tsl'...              % time stamps of detected PSCs that could be fitted
  };

% -- parameters to be computed from raw data and time stamp lists of PSCs
% using algorithms as used in pscfitgui: the prefix 'all' implies being
% determined from all detected (not necessarily fitted) PSCs, whereas
% parameters lacking this suffix are determined from the subset of detected
% AND fitted PSCs
ds.pscDetPar={...
  'allFreq',...         % frequency of detected PSCs (Hz; as determined in threshdetgui)
  'allTsl',...          % time stamps of detected PSCs (ms; as determined in threshdetgui)
  'allAmp',...          % peak amplitude of PSCs (same unit as in raw recording; extracted in pscdeal from raw data in a manner identical to that in PSCFitGui)
  'allTRise20_80',...   % 20-80% (!) rise time (ms; extracted in pscdeal from raw data)
  'tRise20_80',...      % 20-80% (!) rise time of fitted (!) IPSCs (ms; extracted in pscdeal from raw data)
  'thresh'...           % threshold used for detection of PSCs
  };

% -- parameters to be computed from raw data and time stamp lists of PSCs
% via algorithm geared towards dealing with compound PSCs (occurring in
% bursts with at least partly overlapping rise times):
ds.cPscDetPar={...
  'allCAmp',...         % peak amplitude of compound PSCs (same unit as in raw recording)
  'allCTRise',...       % rise time of compound PSCs (ms)  
  'allCTRise20_80',...  % same as above but 20-80%
  };

% -- parameters to be determined from raw data independent of detected or
% fitted PSCs via function phantosic
ds.rawPar={...
  'baseline',...        % baseline (=tonic currents)
  'noise',...           % noise of base line (difference between percentiles, see ds.phPrc)
  'chargePhas'...       % total charge transferred by phasic currents
  };

% -- parameters to be plotted for each experiment
ix=0;
% ix=ix+1;
% ds.plotPar(ix).name='tRise';
% ds.plotPar(ix).bin=[.05:.1:3]; 
% ds.plotPar(ix).bin2=2.^[-3:.4:1]; 
% ix=ix+1;
% ds.plotPar(ix).name='width';
% ds.plotPar(ix).bin=[.25:.25:20]; 
% ds.plotPar(ix).bin2=2.^[0.2:.25:4]; 
% ix=ix+1;
% ds.plotPar(ix).name='allAmp';
% ds.plotPar(ix).bin=[10:5:300 inf];
% ds.plotPar(ix).bin2=[];
% ix=ix+1;
% ds.plotPar(ix).name='allTRise20_80';
% ds.plotPar(ix).bin=[0:.05:5 inf];
% ds.plotPar(ix).bin2=[];

ix=ix+1;
ds.plotPar(ix).name='allCTRise';
ds.plotPar(ix).bin=[0:.05:5 inf];
ds.plotPar(ix).bin2=2.^[-3:.4:1];
ix=ix+1;
ds.plotPar(ix).name='allCAmp';
ds.plotPar(ix).bin=[10:5:300 inf];
ds.plotPar(ix).bin2=2.^(3:.5:7)


% --- settings for raw data analyses
% target sampling freq (Hz)
ds.sampFreq=10000;
% corner frequency of lowpass filter (set to [] for no filtering)
ds.loCFreq=2000;
% settings for differentiator filter (see differfi.m, last three input
% args); set to [] if settings made in threshdetgui shall be applied
% - filter order
ds.differfi.fo=30;
% - lower end of passband (Hz)
ds.differfi.passbf=40;
% - lower end of stopband (Hz, must be substantially higher than
% ds.loCFreq)
ds.differfi.stopbf=3000;
% determine whether computation of PSC amplitudes shall be visualized in a
% separate, temporary figure with interactive features, and upon which
% condition the remaining computations shall be resumed:
% 'none' - no visuals
% 'pause' - computations will resume after figure has been displayed for 1 second
% 'waitfor' - computations will resume once figure has been deleted by user
gr.visualAmpDet='waitfor';

% --- analysis of phasic and tonic currents independent of PSC detection
% and fitting (input parameters into function phantosic):
% - set to positive integer if phantosic computations are to be visualized
% (1=each frame, 2=every second frame, ...; 0 for no visuals)
gr.doMonitorPhantosic=0;
% - segment length (ms)
ds.phIntv=5000;
% - operational method ('peak' or 'Gauss')
ds.phMethod='Gauss';
% - bin width (or bin borders) for amplitude histogram
ds.phBin=2;
% - polarity ('pos' or 'neg'
ds.phPolarity='neg';
% - percentile
ds.phPrc=15.866;

% --- plotting & printing
% set to true if plots of individual experiments are to be produced
gr.doPlot=true;
% graphics format within which to save plot (the usual input arguments into 
% matlab print function; set to [] for no plot)
% gr.printas='-dpsc2';
gr.printas='-djpeg95';
gr.printas=[];
% directory in which to save figures and collected results
gr.fDir=plotPath;

% layout of figure (usual matlab choices 'tall', 'landscape', etc.)
gr.ornt='portrait';
% width of line plot representing data trace
gr.lineW=.5;
% font size (scale bar, title)
gr.fontSz=8;
% scaling factor for post processing, see help of labelscale.m
gr.scaleFac=1;
% marker size
gr.markSz=4;


%% -- 2011_10_13 -------------------------------
expDate='2011_10_13'; 
ds.dDir=[dataPath expDate '_SET3\'];
% file list; column order:
% - file names WITHOUT EXTENSION 
% - independent parameter (e.g. drug concentration)
% - unused
% - unused 
ds.fList={...
  [expDate '_0000'],  0, [], nan;...
  [expDate '_0008'],  1, [], nan;...
  [expDate '_0016'],  2, [], nan;...  
  [expDate '_0024'], 3, [], nan;...
  };
% channel; column order:
% - name
% - unused
ds.chList={'IN 0',[]};

% call pscdeal
pscdeal(ds,gr);

%% -- 2011_10_13 -------------------------------
expDate='2011_10_13'; 
ds.dDir=[dataPath expDate '_SET3\'];
% file list; column order:
% - file names WITHOUT EXTENSION 
% - independent parameter (e.g. drug concentration)
% - unused
% - unused 
ds.fList={...
  [expDate '_0051'],  0, [], nan;...
  [expDate '_0060'],  1, [], nan;...
  [expDate '_0069'],  2, [], nan;...  
  [expDate '_0078'], 3, [], nan;...
  };
% channel; column order:
% - name
% - unused
ds.chList={'IN 0',[]};

% call pscdeal
pscdeal(ds,gr);

%% -- 2011_10_16 -------------------------------
expDate='2011_10_16'; 
ds.dDir=[dataPath expDate '_SET3\'];
% file list; column order:
% - file names WITHOUT EXTENSION 
% - independent parameter (e.g. drug concentration)
% - unused
% - unused 
ds.fList={...
  [expDate '_0000'],  0, [], nan;...
  [expDate '_0008'],  1, [], nan;...
  [expDate '_0016'],  2, [], nan;...  
  [expDate '_0019'], 3, [], nan;...
  };
% channel; column order:
% - name
% - unused
ds.chList={'IN 0',[]};

% call pscdeal
pscdeal(ds,gr);


%% save results
% - retrieve results from figure and save
masterFh=findobj('name','pscdeal master figure','type','figure');
r=get(masterFh,'userdata');
% for compatibility reasons with older versions of pscdeal additionally
% save results fields as variables in capital letters
LISTEXP=r.listExp;
PSCR=r.pscr;
PSCRMN=r.pscrMn;
HISTMONSTER=histogram;
depPar=r.depPar;
save([gr.fDir resultFn],'r','PSCR','PSCRMN','HISTMONSTER','LISTEXP','depPar');
% close master figure 
close(masterFh);

%% simple summary plots
% load data should they have been cleared in the interim
if ~exist('r','var')
  load([gr.fDir resultFn],'r');
end
nSp=numel(r.depPar);
nCol=3;
nRow=ceil(nSp/nCol);

fh=findobj('name',resultFn);
if isempty(fh) 
  fh=mkfig([],'v');
  set(gcf,'name',resultFn,'defaultaxesfontsize',8)
else
  figure(fh)
  clf
end

for g=1:nSp
  subplot(nRow,nCol,g)
  if isfinite(ds.normIpVal)
    % plot of normalized values
    avplot(r.pscrMn(:,:,g),'x',ds.indepPar,'normRow',find(ds.indepPar==ds.normIpVal),...
      'avType','md');
  else
    % plot of absolute values
    avplot(r.pscrMn(:,:,g),'x',ds.indepPar,'avType','md');
    nicexyax(10)
  end
  set(gca,'xticklabel',ds.indepParLabel,'xticklabelrotation',45);
  set(gca,'ygrid','on');
  title(r.depPar{g},'interpreter','none')
end

if ~isempty(gr.printas)
  print(gr.printas,[gr.fDir resultFn '_sumResults']);
end
