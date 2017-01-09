% this is a template of the input structure into pscdeal.m
% -------------------- PREPARATIONS ---------------------------------------
compName=lower(getenv('computername'));
switch compName
  case {'hh-i7'}
    dataP='d:\_data\otc_ctx\ACh\AChBlockIPSC\Ctrl-ACh-Block\';
    plotP='d:\_data\otc_ctx\ACh\AChBlockIPSC\Ctrl-ACh-Block\Figs\';
  case {'hh64','hh-i5'}
    dataP='e:\_data\otc_ctx\ACh\AChBlockIPSC\Ctrl-ACh-Block\';
    plotP='e:\_data\otc_ctx\ACh\AChBlockIPSC\Ctrl-ACh-Block\Figs\';
  case {'eval-lmb'}
    dataP='h:\_data\otc_ctx\ACh\AChBlockIPSC\Ctrl-ACh-Block\';
    plotP='h:\_data\otc_ctx\ACh\AChBlockIPSC\Ctrl-ACh-Block\Figs\';
  otherwise
    error('machine not defined');
end

% --- data set descriptors
% the complete set of values of independent parameter (=concentrations)
ds.indepPar=[0 1 2 3];
% matching descriptors
ds.indepParLabel={'ACh','ACh+Zol','ACh+Dia','ACh (wash)'};

% --- analysis parameters
% -- value of independent parameter corresponding to control (against which
% others will be normalized; set to nan for no normalization)
ds.normIpVal=nan;
% -- value of independent parameter of condition against which others will
% be (statistically) COMPARED
ds.compIpVal=1;

% -- data resulting from pscfitgui
% all fields of struct fitResult to be analyzed (take a look at variable
% fitResult to obtain a complete list)
ds.pscFitPar={'tRise','tDecay','width','amp'};

% -- parameters to be plotted
ds.plotPar(1).name='tRise';
ds.plotPar(1).bin=[.05:.1:3]; 
ds.plotPar(1).bin2=2.^[-3:.4:1]; 
ds.plotPar(2).name='width';
ds.plotPar(2).bin=[.25:.25:14]; 
ds.plotPar(2).bin2=2.^[0.2:.25:4]; 
ds.plotPar(3).name='amp';
ds.plotPar(3).bin=[10:5:300 inf];
ds.plotPar(3).bin2=2.^[4:0.3:8];
ds.plotPar(4).name='allAmp';
ds.plotPar(4).bin=[10:5:300 inf];
ds.plotPar(4).bin2=[];
ds.plotPar(5).name='allTRise20_80';
ds.plotPar(5).bin=[0:.05:8 inf];
ds.plotPar(5).bin2=[];

% -- raw data
% target sampling freq (Hz)
ds.sampFreq=10000;
% interval for xx (ms) (not yet used)
ds.xxIntv=[nan nan];
% if true, analysis of phasic and tonic currents independent of PSC
% detection will be run (via phantosic.m)
ds.phDo=false;
% ** all following parameters are input parameters into function
% phantosic.m:
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

% --- plot appearance 
% width of line plot representing data trace
gr.lineW=.5;
% font size (scale bar, title)
gr.fontSz=8;
% scaling factor for post processing, see help of labelscale.m
gr.scaleFac=1;
% marker size
gr.markSz=4;

% --- layout of figure and plot
% layout of figure (usual matlab choices 'tall', 'landscape', etc.)
gr.ornt='tall';

% --- post-plot action settings
% graphics format within which to save plot (the usual input arguments into 
% matlab print function; set to [] for no plot)
% gr.printas='-dpsc2';
gr.printas='-djpeg95';
gr.printas=[];
% directory to save figures in
gr.fDir=[plotP '\hh\projects\ctx_ACh\'];

% *** global variables that will accumulate results:
% - experiments in rows
% - independent pars in columns
% - parameters in slices
clear global PSCR PSCRMN PSCRVAR LISTEXP HISTMONSTER
global PSCR PSCRMN PSCRVAR LISTEXP HISTMONSTER
PSCR=cell([0 numel(ds.indepPar) numel(ds.pscFitPar)]);
PSCRMN=nan([0 numel(ds.indepPar) numel(ds.pscFitPar)]);
PSCRVAR=nan([0 numel(ds.indepPar) numel(ds.pscFitPar)]);
HISTMONSTER=[];


%% -- 2011_10_13 -------------------------------
expDate='2011_10_13'; 
ds.dDir=[dataP expDate '_SET3\'];
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
depPar=pscdeal(ds,gr);

%% -- 2011_10_13 -------------------------------
expDate='2011_10_13'; 
ds.dDir=[dataP expDate '_SET3\'];
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
depPar=pscdeal(ds,gr);

%% -- 2011_10_16 -------------------------------
expDate='2011_10_16'; 
ds.dDir=[dataP expDate '_SET3\'];
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
depPar=pscdeal(ds,gr);


%% AFTERMATH: quick-and dirty overview plots of averages & variability
figName='AChBlock';
save([gr.fDir figName]);

nSp=numel(depPar);
nCol=3;
nRow=ceil(nSp/nCol);

for h=1:1
  if h==1
    D=PSCRMN;
    ext='_mn';
  else
    D=PSCRVAR./PSCRMN;
    ext='_var';
  end
  figure(100+h), clf, orient landscape
  labelscale('fontSz',8,'scaleFac',1,'lineW',1,'markSz',5);
  for g=1:nSp
    subplot(nRow,nCol,g)
    if isfinite(ds.normIpVal)
      % plot of normalized values
      avplot(D(:,:,g),'x',ds.indepPar,'normRow',find(ds.indepPar==ds.normIpVal),...
        'avType','md');
    else
      % plot of absolute values
      avplot(D(:,:,g),'x',ds.indepPar,'avType','md');
      nicexyax(10)
    end
    set(gca,'xticklabel',ds.indepParLabel);
    set(gca,'ygrid','on');
    title(depPar{g})
  end
  
  if ~isempty(gr.printas),
    print(gr.printas,'-r400',[gr.fDir figName ext]);
  end
end