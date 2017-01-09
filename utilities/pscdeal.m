function depPar=pscdeal(ds_in,gr_in)
% ** function depPar=pscdeal(ds_in,gr_in)
% a utility collecting and processing data produced by PSCFitgui as well as
% (optionally) computing additional parameters from raw data files

% *** HERE'S HOW IT WORKS:
% pscdeal takes two struct arrays. One of them, gr_in, contains graphics
% settings. The other, ds_in, describes analysis parameters and file
% specifics for a set of files representing ONE experiment. The results of
% the data collection and computation are then placed into or appended to
% 3D global variables. By repetitively calling pscdeal, each time with a
% different data set, data can be accumulated.

% *************************************************************************
% This function is a member of a family of functions which perform limited,
% specialized tasks on whole data sets (experiments, each of them composed
% of several files). In contrast to e.g. tslbatch and perievdeal, which are
% called once to work on a whole list of experiments, the functions are
% called for each experiment and place the results (if any) in global
% variables for further use, with varying degrees of safeguards against
% potential problems (duplicate calls on same data, etc.). See the specific
% functions. For orientation, the rough order of evolution is
% [gfplotta,swplotta] -> periburc -> peribuexcit -> pscdeal -> stimdeal
% *************************************************************************

% to do:
% - general overhaul: get rid of PSCRVAR, rethink necessity of global
% variables (persistent var? store in userdata?), computation of 2D
% histograms 
% - the computation of amplitude and rise time of all detected IPSCs (i.e.
% the 'raw' IPSC amplitude, parameters termed 'allAmp' and 'allTRise') is
% lumped with phantosic computations because like phantosic.m it requires
% access to raw data. This does not appear very logical to the end user. In
% a general overhaul, decouple it and make it an 'other par derived from
% raw data';
% - get all parameters into ds (?)
% - set up ds.doPlot

global PSCR PSCRMN PSCRVAR LISTEXP HISTMONSTER
% -------------------------------------------------------------------------
% ------------ PART 0: define defaults
% -------------------------------------------------------------------------
% --- data set descriptors
% the complete set of values of independent parameter (=concentrations)
ds.indepPar=[];
% matching labels
ds.indepParLabel={[]};
% --- analysis parameters
% -- value of independent parameter corresponding to control (against which
% others will be normalized; set to nan for no normalization)
ds.normIpVal=nan;
% -- value of independent parameter of condition against which others will
% be (statistically) COMPARED
ds.compIpVal=[];
% -- data resulting from pscfitgui
% all fields of struct fitResult to be analyzed
ds.pscFitPar={''};
% -- parameters to be plotted (struct)
ds.plotPar=[];
% -- raw data
% target sampling freq (Hz)
ds.sampFreq=nan;
% interval for xx (ms) (not yet used)
ds.xxIntv=[nan nan];
% if true, analysis of phasic and tonic currents independent of PSC
% detection will be run (via phantosic.m)
ds.phDo=false;
% ** all following parameters are input parameters into function
% phantosic.m:
% - segment length (ms)
ds.phIntv=nan;
% - operational method ('peak' or 'Gauss')
ds.phMethod='Gauss';
% - bin width (or bin borders) for amplitude histogram
ds.phBin=nan;
% - polarity ('pos' or 'neg'
ds.phPolarity='neg';
% - percentile
ds.phPrc=nan;

% --- plot appearance 
% width of line plot representing data trace
gr.lineW=.25;
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
gr.doPlot=false;
gr.printas=[];
% directory to save figures in
gr.fDir=[];

% --- file specifics
% data directory
ds.dDir=[''];
% file list; column order:
% - file names WITHOUT EXTENSION 
% - independent parameter (e.g. drug concentration)
% - unused 
% - unused 
ds.fList={...
  'file name', nan, [], [];...
  };
% channel; column order:
% - name
% - unused
ds.chList={'channel name',[]};

% -------------------------------------------------------------------------
% ------------ PART 1: check input, adjust parameters, preparations
% -------------------------------------------------------------------------

ds=checkFields(ds,ds_in);
gr=checkFields(gr,gr_in);  

% --- graphics
labelscale('fontSz',gr.fontSz,'scaleFac',gr.scaleFac,'lineW',gr.lineW,'markSz',gr.markSz); 
nFile=size(ds.fList,1);

% main figure
switch gr.ornt
  case 'landscape'
    fh1=mkfig(1,'b');
  case 'tall'
    fh1=mkfig(1,'v');
end
clf
orient(gr.ornt);

% colormap(bone)

% figure name: first of file names in set & channel name(s)
tmp=[ds.chList{:,1}];
figName=[ds.fList{1,1} '_' tmp(~isspace(tmp))];

% number of indpendent and dependent parameters
nIndepPar=numel(ds.indepPar);
nPscFitPar=numel(ds.pscFitPar);
% additional parameters created in the code further below:
% - time stamps of detected PSCs
% - time stamps of those that could be fitted
% - frequency of detected PSCs
% - frequency of those that could be fitted
% - charge transferred per psc (of those that could be fitted)
% - total charge transferred (by those that could be fitted)
otherPar={'allTsl','tsl','freq','freqFit','chargePPsc','chargePscTot'};
nOtherPar=numel(otherPar);
if ds.phDo
  % parameters extracted from raw data, mostly produced by phantosic: 
  % - base line (=tonic currents)
  % - noise of base line
  % - total charge transferred by phasic currents
  % - amplitude of all detected (not necessarily fitted) IPSCs (outside phantosic)
  % - 20-80 % (!) rise time of all detected (not necessarily fitted) IPSCs (outside phantosic)
  % - 20-80 % (!) rise time of fitted IPSCs (outside phantosic)  
  phPar={'baseline','noise','chargePhas','allAmp','allTRise20_80','tRise20_80'};
else
  phPar={};
end
nPhPar=numel(phPar);
% so, total set of dependent parameters: the ones listed in ds.pscFitPar
% first, then derived psc pars, then phantosic pars 
depPar=cat(2,ds.pscFitPar,otherPar,phPar);
nDepPar=numel(depPar);

% check ds.indepPar
if ~isempty(setdiff([ds.fList{:,2}],ds.indepPar))
  error('check independent parameters');
end

% initialize LISTEXP
if ~istable(LISTEXP)
  tmpCell=cell(1,numel(ds.indepPar));
  [tmpCell{:}]=deal('');
  LISTEXP=cell2table(tmpCell);
  LISTEXP.Properties.VariableNames=strrep(strrep(ds.indepParLabel,'+','p'),'-','m');
end  

% 2D histogram template
histTemplate=nan(numel(ds.plotPar(1).bin2),numel(ds.plotPar(2).bin2));

% -------------------------------------------------------------------------
% ------------ PART 2: loop over files & load and process data
% -------------------------------------------------------------------------

for g=1:nFile
  % clear anything that may have remained from previous file
  clear global bu evt 
  clear fitResult fitHead head
  % index into column of global results variable PSCR
  colIx=find(ds.indepPar==ds.fList{g,2});
  % channel name
  if g==1
    deblChName=ds.chList{:,1};
    deblChName=deblChName(~isspace(deblChName));
    % name of current experiment: first file + channel name
    curExpName=[ds.fList{1,1} '_' deblChName];
    % row index for current experiment into global results variable: if
    % experiment name exists already, point to corresponding row
    rowIx=find(strcmp(curExpName,LISTEXP.Properties.RowNames));
    %     rowIx=find(strcmp(curExpName,LISTEXP));
    if isempty(rowIx)
      rowIx=size(PSCR,1)+1;
      [~,ia]=intersect(ds.indepPar,cat(2,ds.fList{:,2}),'stable');
      LISTEXP{rowIx,ia}=ds.fList(:,1)';
      LISTEXP.Properties.RowNames{rowIx}=curExpName;
      % ** preallocate
      PSCRMN(rowIx,:,1:nDepPar)=nan;
      PSCRVAR(rowIx,:,1:nDepPar)=nan;
      [PSCR(rowIx,:,1:nDepPar)]=deal({[nan]});
      % dimensions: first par, second par, experiment, indep par value
      HISTMONSTER(1:numel(ds.plotPar(1).bin2),1:numel(ds.plotPar(2).bin2),...
        rowIx,1:nIndepPar)=repmat(histTemplate,[1 1 1 nIndepPar]);
    elseif numel(rowIx)>1
      error('internal:duplicate experiments in LISTEXP');
    else
      disp(['overwriting values for ' curExpName]);
      % ** wipe any existing entries 
      PSCRMN(rowIx,:,1:nDepPar)=nan;
      PSCRVAR(rowIx,:,1:nDepPar)=nan;
      [PSCR(rowIx,:,1:nDepPar)]=deal({[nan]});
      % ��� [HISTMONSTER(1,:)]=deal({histTemplate});
    end
  end
  % base name of processed results file
  fnChnBase=[ds.fList{g,1} '_' deblChName];
  fn=ds.fList{g,1};
  pscFn=([ds.dDir fnChnBase '_IPSC_res.mat']);
  
  % ��
  
  pscFn=([ds.dDir fnChnBase '_EPSC_res.mat']);
  
  % load data (most importantly, fit results and original tsl)
  load(pscFn)
  for idi=1:nPscFitPar
    switch(ds.pscFitPar{idi})
      case 'width'
        if ~isfield(fitResult,'width')
          % *** (half-width) has been added to PSCFit in May 2012, i.e. it
          % must be computed for data before this date
          if size(fitResult.tDecay,1)==1
            % if we're dealing with only one component half width is
            % tau*ln(2)
            fitResult.width=fitResult.tDecay(end,:)*log(2);
          else
            nCutout=size(fitResult.OKIx,2);
            fitResult.width=nan(1,nCutout);
            % �� should be made a parameter
            widthFac=.5;
            % discrete time
            ix3=(0:fitHead.wp.xIntvFitEndMax_pts);
            for cIx=1:nCutout
              if fitResult.OKIx(1,cIx)
                a1=fitResult.amp(1,cIx);
                a2=fitResult.amp(2,cIx);
                tau1=fitResult.tDecay(1,cIx)/(fitHead.wp.si/1000);
                tau2=fitResult.tDecay(2,cIx)/(fitHead.wp.si/1000);
                cutoutFit=a1./exp(ix3/tau1)+a2./exp(ix3/tau2);
                % - indexes to points of fit flanking intersection with zero
                x1=find(cutoutFit>(a1+a2)*widthFac,1,'last');
                if ~isempty(x1) && x1<fitHead.wp.xIntvFitEndMax_pts
                  x2=x1+1;
                  % corresponding amplitudes
                  y1=cutoutFit(x1);
                  y2=cutoutFit(x2);
                  % intersection:
                  fitResult.width(1,cIx)=x1+((a1+a2)*widthFac-y1)*(x2-x1)/(y2-y1);
                end
              end
            end
            fitResult.width=fitResult.width*fitHead.wp.si/1000;
          end
        end
        % proceed as usual:
        % copy individual IPSC data into corresponding cell
        PSCR{rowIx,colIx,idi}=fitResult.(ds.pscFitPar{idi})(end,:)';
        % copy average (=median!) of these into array
        PSCRMN(rowIx,colIx,idi)=nanmedian(fitResult.(ds.pscFitPar{idi})(end,:)');
        PSCRVAR(rowIx,colIx,idi)=diff(prctile(fitResult.(ds.pscFitPar{idi})(end,:)',[2.5 97.5]));
      case 'tDecay'
        % pick fast, slow or weighted:
        % copy individual IPSC data into corresponding cell
        PSCR{rowIx,colIx,idi}=fitResult.(ds.pscFitPar{idi})(end,:)';
        % copy average (=median!) of these into array
        PSCRMN(rowIx,colIx,idi)=nanmedian(fitResult.(ds.pscFitPar{idi})(end,:)');
        PSCRVAR(rowIx,colIx,idi)=diff(prctile(fitResult.(ds.pscFitPar{idi})(end,:)',[2.5 97.5]));
      otherwise
        % all other cases not requiring special treatment:
        % copy individual IPSC data into corresponding cell
        PSCR{rowIx,colIx,idi}=fitResult.(ds.pscFitPar{idi})(end,:)';
        % copy average (=median!) of these into array
        PSCRMN(rowIx,colIx,idi)=nanmedian(fitResult.(ds.pscFitPar{idi})(end,:)');
        PSCRVAR(rowIx,colIx,idi)=diff(prctile(fitResult.(ds.pscFitPar{idi})(end,:)',[2.5 97.5]));
    end
  end

  % additional parameters
  for idi=1:nOtherPar
    offs=nPscFitPar;
    % variability not defined or not requested for most of the parameters
    % dealt with here
    PSCRVAR(rowIx,colIx,idi+offs)=nan;
    switch otherPar{idi}
      case 'tsl'
        if exist('fitResult','var')
          % ts of PSCs which were fitted, which must be a subset of all ts -
          % check this
          [fitTsl,fitTsIx]=intersect(evt.tsl{1},fitResult.tsl);
          if ~isequal(fitTsl(:),fitResult.tsl(:))
            % see whether they really differ
            if any(abs(fitTsl(:)-fitResult.tsl(:))>.1)
              error('tsl of fitted events is not a complete subset of all detected events');
            end
          end
          PSCR{rowIx,colIx,idi+offs}=fitResult.tsl;
        end
        PSCRMN(rowIx,colIx,idi+offs)=nan;
        
      case 'allTsl'
        % simply stuff tsl in there
        PSCR{rowIx,colIx,idi+offs}=evt.tsl{1};
        PSCRMN(rowIx,colIx,idi+offs)=nan;
      case 'freq'
        % frequency of detected (not necessarily fitted) PSCs in Hz
        f=numel(evt.tsl{1})/diff(head.ds.fileInfo.recTime);
        PSCR{rowIx,colIx,idi+offs}=f;
        PSCRMN(rowIx,colIx,idi+offs)=f;
      case 'freqFit'
        % frequency of fitted PSCs in Hz
        % �� this should be superfluous with the version of pscfitgui 2.1
        % as of Oct 30, 2011
        if exist('fitResult','var')
          tsl=fitResult.tsl(isfinite(fitResult.amp(end,:)));
          f=numel(tsl)/diff(head.ds.fileInfo.recTime);
          PSCR{rowIx,colIx,idi+offs}=f;
          PSCRMN(rowIx,colIx,idi+offs)=f;
        end
      case 'chargePPsc'
        % integral/charge transferred per psc as computed on the basis of
        % fitted parameters in original current units/s, so for pA it's
        % 10^-12 C
        if exist('fitResult','var')
          c=fitResult.amp(end,:).*fitResult.tDecay(end,:)/1000;
          PSCR{rowIx,colIx,idi+offs}=c';
          PSCRMN(rowIx,colIx,idi+offs)=nanmedian(c);
          PSCRVAR(rowIx,colIx,idi+offs)=diff(prctile(c,[2.5 97.5]));
        end
      case 'chargePscTot'        
        % total charge transferred by PSCs (as computed on the basis of
        % fitted parameters), divided by recording time (so, this is the
        % average current in the original units)
        if exist('fitResult','var')
          c=nansum(fitResult.amp(end,:).*fitResult.tDecay(end,:)/1000)/diff(head.ds.fileInfo.recTime);
          PSCR{rowIx,colIx,idi+offs}=c;
          PSCRMN(rowIx,colIx,idi+offs)=c;
        end
      otherwise 
        error('bad derived psc par');
    end
  end
  
  %%% computations requiring raw data, including phantosic:
  if ds.phDo
    % raw data
    rawFn=[ds.dDir '\' fn '.abf'];
    if exist(rawFn,'file')
      [d,si,fi]=abfload(rawFn,'channels',ds.chList(:,1)');
      % lowpass filter
      if isfield(fitHead.ap,'loCFreq') && isfinite(fitHead.ap.loCFreq)
        d=lofi(d,si,fitHead.ap.loCFreq);
      end
      % --- phantosic 
      nPts=ds.phIntv/(si/1e3);
      [base,dev,phas,gof]=phantosic(d,nPts,ds.phBin,'method',ds.phMethod,...
        'polarity',ds.phPolarity,'prc',ds.phPrc,'frame',0,'pau',1);
      
      % --- peak and rise time
      % this is a version of the peak detection algorithm derived from the
      % one in pscfitguifunc and mostly identical with it
      d=sgolayfilt(d,fitHead.ap.smoothPolyOrder,fitHead.wp.smoothTSpan_pts);
      % generate cutouts: beginning of base line interval to end of
      % peak det interval
      [cutout,isCutout]=tsl2exc(d,si,evt.tsl,'win',[fitHead.ap.xIntvBaseLine(1) fitHead.ap.xIntvPeak(2)]);
      % 'no event left behind (?)': as the current cutout window may differ
      % from the original one make sure that lost cutouts are also removed
      % from etsl and that - further below - they are properly embedded in
      % an array (in PSCR) of the same size as evt.tsl
      tsl=evt.tsl{1}(isCutout{1});
      % size of things
      [coNRow,coNCol]=size(cutout);
      % convert time intervals to pts:
      % - base line
      baseIx=cont2discrete([0 diff(fitHead.ap.xIntvBaseLine)],si/1000,'intv',1);
      % - peak search interval
      peakIx=cont2discrete(fitHead.ap.xIntvPeak-fitHead.ap.xIntvBaseLine(1),si/1000,'intv',1);
      % definition of interval to look for rise time as in pscfitguifunc
      riseTIx=cont2discrete([fitHead.ap.xIntvBaseLine(2) fitHead.ap.xIntvPeak(2)]-fitHead.ap.xIntvBaseLine(1),si/1000,'intv',1);
      % subtract base line from traces
      cutout=cutout-repmat(mean(cutout(baseIx(1):baseIx(2),:)),coNRow,1);
      % invert?
      if strcmp(ds.phPolarity,'neg')
        cutout= -cutout;
      end
      % determine peak: 
      % - restrict individual search windows to ocurrence of next detected
      % PSC, shifted a bit to the left (so, here is the index to the last
      % point within cutouts within which to look for peak):
      lastGoodIx=min(peakIx(2),cont2discrete([diff(tsl); fitHead.ap.xIntvPeak(2)]-fitHead.ap.xIntvBaseLine(1)-.5,si/1000));
      % - preallocate
      pk=nan(1,coNCol);
      for evIx=1:coNCol
        % - detect all positive peaks
        excIx=peakIx(1):lastGoodIx(evIx);
        dfC=diff(cutout(excIx,evIx));
        pkIx=find((diff(dfC<=0)==1))+1;
        if ~isempty(pkIx)
          % - find maximal one
          [m,mIx]=max(cutout(excIx(pkIx),evIx));
          % - average of maximal peak and neighboring points 
          pk(evIx)=mean(cutout(excIx(pkIx(mIx))+[-1:1],evIx));
          % plot to check/debug
          if 0
            plot(cutout(:,evIx),'k')
            hold on
            ph=plot(lastGoodIx(evIx),0,'rv');
            ph=plot(excIx(pkIx(mIx)),pk(evIx),'co');
            hold off
          end
        end
      end
      % 20-80% rise time: upsample unfiltered waveforms to 100 kHz (si=10
      % us) for more precision (edge effects of built-in filter of function
      % resample should not play a role as baseline is subtracted already
      % and the right edge, prone to distortion, is far beyond the peak)
      iSi=10;
      iCutout=resample(cutout(riseTIx(1):riseTIx(2),:),si,iSi);
      niSample=size(iCutout,1);
      [~,i1]=max(iCutout>=repmat(.2*pk,niSample,1));
      [~,i2]=max(iCutout>=repmat(.8*pk,niSample,1));
      tRise=(i2-i1)'*(iSi/1000);
      % events with no peak will show up as zeroes in the rise time, so 
      % do something about it
      tRise(~isfinite(pk))=nan;
      
      if 0 && any(tRise<.1)
        figure(1)
        clf
        hold on
        fastIx=tRise<.1;
        plot(cutout(:,fastIx),'o-')
        niceyax;
      end
      
      % if any of the cutouts could not be produced, make sure data are
      % embedded properly
      if any(~isCutout{1})
        disp('at least one cutout could not be reproduced - embedding data');
        templateArr=nan(numel(evt.tsl{1}),1);
        maxPeak=templateArr;
        maxPeak(isCutout{1})=pk';
        % use templateArr as temp var for tRise
        templateArr(isCutout{1})=tRise;
        tRise=templateArr;
      else
        maxPeak=pk';
        % tRise is fine
      end
      
    else
      warning([rawFn ' does not exist']);
    end
  end
  
  
   %%% embed phantosic and other params
   for idi=1:nPhPar
    offs=nPscFitPar+nOtherPar;
    % variability not defined or not requested for the parameters dealt
    % with here
    PSCRVAR(rowIx,colIx,idi+offs)=nan;
    switch phPar{idi}
      case 'baseline'
        % baseline for recording
        PSCR{rowIx,colIx,idi+offs}=base;
        PSCRMN(rowIx,colIx,idi+offs)=nanmean(base);
      case 'noise'
        % deviation from baseline (noise)
        PSCR{rowIx,colIx,idi+offs}=dev;
        PSCRMN(rowIx,colIx,idi+offs)=nanmean(dev);
      case 'chargePhas'
        % phasic current
        PSCR{rowIx,colIx,idi+offs}=phas;
        PSCRMN(rowIx,colIx,idi+offs)=nanmean(phas);
      case 'allAmp'
        % stuff into cell arrays
        PSCR{rowIx,colIx,idi+offs}=maxPeak;
        PSCRMN(rowIx,colIx,idi+offs)=nanmedian(maxPeak);
      case 'allTRise20_80'
        % stuff into cell arrays
        PSCR{rowIx,colIx,idi+offs}=tRise;
        PSCRMN(rowIx,colIx,idi+offs)=nanmedian(tRise);
      case 'tRise20_80'
        % stuff rise time of identified PSCs (see above) into cell arrays
        PSCR{rowIx,colIx,idi+offs}=tRise(fitTsIx);
        PSCRMN(rowIx,colIx,idi+offs)=nanmedian(tRise(fitTsIx));
      otherwise
        error('bad derived psc par');
    end
   end
  
  % computations for 2D histogram
  % - transform PSCR & extract parameters of interest
  tmppscd=permute(PSCR(rowIx,:,:),[1 3 2]);
  pscd=[tmppscd{1,strcmp(ds.plotPar(1).name,[ds.pscFitPar otherPar phPar]),colIx},...
    tmppscd{1,strcmp(ds.plotPar(2).name,[ds.pscFitPar otherPar phPar]),colIx}];
  n=hist3(pscd,'Edges',{ds.plotPar(1).bin2 ds.plotPar(2).bin2});
  % - stuff histogram into appropriate place in HISTMONSTER
  HISTMONSTER(:,:,rowIx,colIx)=n;
   
end


% -------------------------------------------------------------------------
% ------------ PART 3: plot data
% -------------------------------------------------------------------------

% currently, if ds.printas is empty, don't plot
if gr.doPlot
  
  % 1. -------- univariate cumulative histograms of IPSC parameters ---------
  figure(fh1); clf,
  set(gcf,'defaultlinelinewidth',1.5);
  % colors to be used for display of several histograms
  stairsCm=[...
    .9       0      0;...
    0       0      0;...
    0.0667    0.2582    0.8165;...
    0       .8      0;...
    .7        0         .7;...
    0.0100    0.3162    0.3162;...
    ];
  stairsCm=repmat(stairsCm,[ceil(nIndepPar/size(stairsCm,1)) 1]);
  pars={ds.plotPar.name};
  nPars=numel(pars);
  
  for parIx=1:nPars
    hsp=subplot(max(3,nPars),1,parIx);
    hold on
    % ** rowIx has been computed above and can still be used
    for g=1:nFile
      % index into column of global results variable PSCR
      colIx=find(ds.indepPar==ds.fList{g,2});
      % pull out data
      helpIx=strcmp(depPar,pars{parIx});
      n=PSCR{rowIx,colIx,helpIx};
      hh=histc(n,ds.plotPar(parIx).bin);
      % ** normalize to number of events
      hh_norm=hh/sum(hh);
      hs=stairs(ds.plotPar(parIx).bin,hh_norm);
      set(hs,'Color',stairsCm(g,:));
    end
    grid on
    axis([ds.plotPar(parIx).bin([1 end]) -.05 1.05])
    xlabel(ds.plotPar(parIx).name);
    ylabel('N (norm.)');
  end
  
  legend(ds.indepParLabel);
  
  if ~isempty(gr.printas),
    print(gr.printas,'-r400',[gr.fDir figName '_hist']);
  end
  drawnow
  pause(.05);
  
  % 2. ------------ 2D histograms of IPSC parameters -------------
  fh2=mkfig(2,'b'); clf, orient landscape;
  colormap(coma('CubicL'));
  % layout of subplots
  nRow=2;
  nCol=numel(ds.indepPar);
  % local normalization to number of PSCs?
  doNormalize=false;
  
  % ** rowIx=index for current experiment into global results variable; has
  % been computed above and can still be used
  
  % start by pulling out control histogram
  compColIx=find(ds.indepPar==ds.compIpVal);
  % pull out histogram
  compHist=HISTMONSTER(:,:,rowIx,compColIx);
  if doNormalize
    % normalize
    compHist=compHist/sum(compHist(:));
  end
  
  sph=nan(2,nFile);
  for g=1:nFile
    % index into column of global results variable PSCR
    colIx=find(ds.indepPar==ds.fList{g,2});
    % pull out histogram
    hi=HISTMONSTER(:,:,rowIx,colIx);
    if doNormalize
      % normalize
      hi=hi/sum(hi(:));
    end
    cLim(g)=prctile(hi(:),99);
    % plot
    sph(1,g)=subplot(nRow,nCol,g);
    % hi has to be flipped in order for the first plot parameter to be
    % plotted along the abscissa
    [~,h]=contourf(ds.plotPar(1).bin2,ds.plotPar(2).bin2,hi',30);
    set(h,'linestyle','none');
    grid on
    xlabel(ds.plotPar(1).name);
    ylabel(ds.plotPar(2).name);
    title(ds.indepParLabel{g});
    % plot difference
    sph(2,g)=subplot(nRow,nCol,nCol+g);
    axis off
    if g==compColIx+1
      [~,h]=contourf(ds.plotPar(1).bin2,ds.plotPar(2).bin2,hi'-compHist',30);
      axis on
      set(h,'linestyle','none');
      grid on
      xlabel(ds.plotPar(1).name);
      ylabel(ds.plotPar(2).name);
      title([ds.indepParLabel{g} ' - ' ds.indepParLabel{compColIx}]);
    end
  end
  set(sph(1,:),'clim',[0 max(cLim)]);
  set(sph(2,:),'clim',max(cLim)*[-1 1]);
  drawnow
  
  if ~isempty(gr.printas),
    print(gr.printas,'-r400',[gr.fDir figName '_2Dhist_' ds.plotPar(1).name '_' ds.plotPar(2).name]);
  end
  
end
drawnow
pause(.05);


% ------------ LOCAL FUNCTIONS ------------------------------------------
% ------------ LOCAL FUNCTIONS ------------------------------------------
function f=checkFields(f,f_in)
s=fieldnames(f);
s_in=fieldnames(f_in);
sect=setdiff(s_in,s);
if ~isempty(sect)
  errordlg({'nonmatching field names(s): '; strvcat(sect)});
  error('see error dialogue');
end
% assign input values of input struct to current one
for g=1:numel(s_in)
  f.(s_in{g})=f_in.(s_in{g});
end
