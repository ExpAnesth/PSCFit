function ud=pscdeal(ds_in,gr_in)
% ** function ud=pscdeal(ds_in,gr_in)
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
% - update computation of 2D histograms
% - get all parameters into ds (?)

% =========================================================================
%          PART 1: check input, adjust parameters, preparations
% =========================================================================
[ds,gr]=pscdeal_defaultParams;
ds=checkFields(ds,ds_in);
gr=checkFields(gr,gr_in);
nFile=size(ds.fList,1);

% figure name: first of file names in set & channel name(s)
tmp=[ds.chList{:,1}];
figName=[ds.fList{1,1} '_' tmp(~isspace(tmp))];

% number of indpendent and dependent parameters
nIndepPar=numel(ds.indepPar);
nPscFitPar=numel(ds.pscFitPar);
nPscDetPar=numel(ds.pscDetPar);
nCPscDetPar=numel(ds.cPscDetPar);
nRawPar=numel(ds.rawPar);

% so, total set of dependent parameters: the ones listed in ds.pscFitPar
% first, then ds.pscDetPars, then ds.cPscDetPars, then raw pars
depPar=cat(2,ds.pscFitPar,ds.pscDetPar,ds.cPscDetPar,ds.rawPar);
nDepPar=numel(depPar);

% find out whether any of the requested parameters requires upload of raw
% data
if any(ismember(depPar,{'allAmp','allTRise20_80','tRise20_80',...
    'allCAmp','allCTRise','allCTRise20_80',...
    'baseline','noise','chargePhas'}))
  doReadRawData=true;
else
  doReadRawData=false;
end
% 'allAmp','allTRise20_80','tRise20_80' are computed in a manner identical
% to that in pscfitgui, so we need analysis parameters (base line window
% etc.) as used in pscfitgui
if nPscFitPar>0 || any(ismember(depPar,{'allAmp','allTRise20_80','tRise20_80'})) 
  doNeedFitResults=true;
else
  doNeedFitResults=false;
end
% check whether bivariate histogram of PSC parameters shall be computed
if numel(ds.plotPar)>1
  doBivarHistogram=true;
else
  doBivarHistogram=false;
end
% check ds.indepPar
if ~isempty(setdiff([ds.fList{:,2}],ds.indepPar))
  error('check independent parameters');
end

% bivariate histogram template
if doBivarHistogram
  histTemplate=nan(numel(ds.plotPar(1).bin2),numel(ds.plotPar(2).bin2));
end

switch gr.visualAmpDet
  case 'none' 
    doMonitorAmpDet=false;
  case {'pause','waitfor'}
    doMonitorAmpDet=true;
  otherwise
    error('bad value for gr.visualAmpDet')
end

% find/set up master figure
masterFh=findobj('name','pscdeal master figure','type','figure');
if isempty(masterFh)
  masterFh=mkfig(1,'min');
  set(masterFh,'name','pscdeal master figure','units','normalized',...
    'toolbar','none','menubar','none');
  uicontrol('Units','normalized','position',[.1 .8 .5 .1],...
    'style','text','string','Analyzed experiments:');
  listH=uicontrol('Units','normalized','position',[.1 .1 .5 .7],...
    'style','listbox','tag','listexp');
  % set up userdata
  ud.depPar=depPar;
  % - for compatibility with pre-table Matlab and older versions of pscdeal
  % set up the list of experiment names and matching files as separate cell
  % arrays
  ud.listExp=cell(1);
  ud.listExpFile=cell(1,numel(ds.indepPar));
  ud.pscr=cell([0 numel(ds.indepPar) numel(ds.pscFitPar)]);
  ud.pscrMn=nan([0 numel(ds.indepPar) numel(ds.pscFitPar)]);
  ud.histogram=[];
else
  ud=get(masterFh,'userdata');
  listH=findall(masterFh,'tag','listexp');
end

% =========================================================================
%            PART 2: loop over files & load and process data
% =========================================================================

for g=1:nFile
  disp(' ');
  disp('------------------------------------------------------------------------------')
  disp(['processing ' ds.fList{g,1} '...']);
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
    % row index for current experiment into results variables: if
    % experiment name exists already, point to corresponding row
    rowIx=find(strcmp(curExpName,ud.listExp));
    
    if isempty(rowIx)
      rowIx=size(ud.pscr,1)+1;
      [~,ia]=intersect(ds.indepPar,cat(2,ds.fList{:,2}),'stable');
      ud.listExpFile(rowIx,ia)=ds.fList(:,1)';
      ud.listExp{rowIx,1}=curExpName;
      % ** preallocate
      ud.pscrMn(rowIx,:,1:nDepPar)=nan;
      [ud.pscr(rowIx,:,1:nDepPar)]=deal({nan});
      if doBivarHistogram
        % histogram dimensions: first par, second par, experiment, indep par value
        ud.histogram(1:numel(ds.plotPar(1).bin2),1:numel(ds.plotPar(2).bin2),...
          rowIx,1:nIndepPar)=repmat(histTemplate,[1 1 1 nIndepPar]);
      end
    elseif numel(rowIx)>1
      error('internal:duplicate experiments in ud.listExp');
    else
      disp('** overwriting values');
      % ** wipe any existing entries
      ud.pscrMn(rowIx,:,1:nDepPar)=nan;
      [ud.pscr(rowIx,:,1:nDepPar)]=deal({nan});
      % §§§ [ud.histogram(1,:)]=deal({histTemplate});
    end
  end
  % base name of processed results file
  fnChnBase=[ds.fList{g,1} '_' deblChName];
  fn=ds.fList{g,1};
  pscFn=([ds.dDir fnChnBase ds.fnSuffix '.mat']);
  % load data (most importantly, fit results and original tsl)
  load(pscFn)
  existFitResult=logical(exist('fitResult','var'));
  if doNeedFitResults && ~existFitResult
    error('Fit parameters are requested, but current data were not analyzed by PSCFit');
  end
  % -----------------------------------------------------------------------  
  %             collect/compute parameters of fitted PSCs 
  % -----------------------------------------------------------------------
  disp('collecting & computing parameters of fitted PSCs...')
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
            % §§ should be made a parameter
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
        ud.pscr{rowIx,colIx,idi}=fitResult.(ds.pscFitPar{idi})(end,:)';
        % copy average (=median!) of these into array
        ud.pscrMn(rowIx,colIx,idi)=nanmedian(fitResult.(ds.pscFitPar{idi})(end,:)');
      case 'tDecay'
        % pick fast, slow or weighted:
        % copy individual IPSC data into corresponding cell
        ud.pscr{rowIx,colIx,idi}=fitResult.(ds.pscFitPar{idi})(end,:)';
        % copy average (=median!) of these into array
        ud.pscrMn(rowIx,colIx,idi)=nanmedian(fitResult.(ds.pscFitPar{idi})(end,:)');
      case 'qFit'
        c=fitResult.qFit;
        ud.pscr{rowIx,colIx,idi}=c';
        ud.pscrMn(rowIx,colIx,idi)=nanmedian(c);
      case 'tsl'
        % ts of PSCs which were fitted, which must be a subset of all ts -
        % check this
        [fitTsl,fitTsIx]=intersect(evt.tsl{1},fitResult.tsl);
        if ~isequal(fitTsl(:),fitResult.tsl(:))
          % see whether they really differ
          if any(abs(fitTsl(:)-fitResult.tsl(:))>.1)
            error('tsl of fitted events is not a complete subset of all detected events');
          end
        end
        ud.pscr{rowIx,colIx,idi}=fitResult.tsl(:);
        ud.pscrMn(rowIx,colIx,idi)=nan;
      case 'freqFit'
        % frequency of fitted PSCs in Hz
        % §§ this should be superfluous with the version of pscfitgui 2.1
        % as of Oct 30, 2011
        tsl=fitResult.tsl(isfinite(fitResult.amp(end,:)));
        f=numel(tsl)/diff(head.ds.fileInfo.recTime);
        ud.pscr{rowIx,colIx,idi}=f;
        ud.pscrMn(rowIx,colIx,idi)=f;
      case 'chargePPsc'
        % integral/charge transferred per psc as computed on the basis of
        % fitted parameters in original current units/s, so for pA it's
        % 10^-12 C
        c=fitResult.amp(end,:).*fitResult.tDecay(end,:)/1000;
        ud.pscr{rowIx,colIx,idi}=c';
        ud.pscrMn(rowIx,colIx,idi)=nanmedian(c);
      case 'chargePscTot'
        % total charge transferred by PSCs (as computed on the basis of
        % fitted parameters), divided by recording time (so, this is the
        % average current in the original units)
        c=nansum(fitResult.amp(end,:).*fitResult.tDecay(end,:)/1000)/diff(head.ds.fileInfo.recTime);
        ud.pscr{rowIx,colIx,idi}=c;
        ud.pscrMn(rowIx,colIx,idi)=c;
      otherwise
        % all other cases not requiring special treatment:
        % copy individual IPSC data into corresponding cell
        ud.pscr{rowIx,colIx,idi}=fitResult.(ds.pscFitPar{idi})(end,:)';
        % copy average (=median!) of these into array
        ud.pscrMn(rowIx,colIx,idi)=nanmedian(fitResult.(ds.pscFitPar{idi})(end,:)');
    end
  end
  
  % -----------------------------------------------------------------------
  %                  read raw data if required 
  % -----------------------------------------------------------------------
  if doReadRawData
    % raw data
    rawFn=[ds.dDir '\' fn '.abf'];
    if exist(rawFn,'file')
      [d,si]=abfload(rawFn,'channels',ds.chList(:,1)');
    else
      warning([rawFn ' does not exist']);
    end
  end
  
  % -----------------------------------------------------------------------
  % collect parameters of detected PSCs and compute amplitudes and rise
  % times of from raw data, pscfitgui-style
  % -----------------------------------------------------------------------
  % loCFreq is the corner frequency of the lowpass which d is run through
  % to yield loD, which may be used further below
  loCFreq=[];
  if any(ismember(depPar,{'allAmp','allTRise20_80','tRise20_80'}))
    disp('computing amplitudes and tRise of detected PSCs (pscfitgui algorithm)...')
    % set filter frequency as used in pscfitgui, otherwise set 'global'
    % freq
    if isfield(fitHead.ap,'loCFreq') && isfinite(fitHead.ap.loCFreq)
      loCFreq=fitHead.ap.loCFreq;
    else
      loCFreq=ds.loCFreq;
    end
    if ~isempty(loCFreq)
      loD=lofi(d,si,loCFreq);
    else
      loD=d;
    end
    smoothD=sgolayfilt(loD,fitHead.ap.smoothPolyOrder,fitHead.wp.smoothTSpan_pts);
    % generate cutouts: beginning of base line interval to end of
    % peak det interval
    [cutout,isCutout]=tsl2exc(smoothD,si,evt.tsl,'win',[fitHead.ap.xIntvBaseLine(1) fitHead.ap.xIntvPeak(2)]);
    % 'no event left behind (?)': as the current cutout window may differ
    % from the original one make sure that lost cutouts are also removed
    % from etsl and that - further below - they are properly embedded in
    % an array (in ud.pscr) of the same size as evt.tsl
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
        [~,mIx]=max(cutout(excIx(pkIx),evIx));
        % - average of maximal peak and neighboring points
        pk(evIx)=mean(cutout(excIx(pkIx(mIx))+(-1:1),evIx));
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
      figure(10)
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
  end
  % embed parameters, some of them possibly computed above
  for idi=1:nPscDetPar
    sliceIx=strcmp(ds.pscDetPar{idi},depPar);
    switch ds.pscDetPar{idi}
      case 'allFreq'
        % frequency of detected (not necessarily fitted) PSCs in Hz
        f=numel(evt.tsl{1})/diff(head.ds.fileInfo.recTime);
        ud.pscr{rowIx,colIx,sliceIx}=f;
        ud.pscrMn(rowIx,colIx,sliceIx)=f;
      case 'allTsl'
        ud.pscr{rowIx,colIx,sliceIx}=evt.tsl{1};
        ud.pscrMn(rowIx,colIx,sliceIx)=nan;
      case 'allAmp'
        ud.pscr{rowIx,colIx,sliceIx}=maxPeak;
        ud.pscrMn(rowIx,colIx,sliceIx)=nanmedian(maxPeak);
      case 'allTRise20_80'
        ud.pscr{rowIx,colIx,sliceIx}=tRise;
        ud.pscrMn(rowIx,colIx,sliceIx)=nanmedian(tRise);
      case 'tRise20_80'
        % rise time of identified PSCs (see above) 
        ud.pscr{rowIx,colIx,sliceIx}=tRise(fitTsIx);
        ud.pscrMn(rowIx,colIx,sliceIx)=nanmedian(tRise(fitTsIx));
      case 'thresh'
        ud.pscr{rowIx,colIx,sliceIx}=head.ap.thresh;
        ud.pscrMn(rowIx,colIx,sliceIx)=head.ap.thresh;
      otherwise
        error('bad pscDet par');
    end
  end
  
  % -----------------------------------------------------------------------
  % compute amplitudes and rise times of detected PSCs from raw data via
  % algorithm designed for compound PSCs
  % -----------------------------------------------------------------------  
  if any(ismember(depPar,{'allCAmp','allCTRise','allCTRise20_80'}))
    disp('computing amplitudes and tRise of detected PSCs (''compound'' PSC algorithm)...')
    cAmpDetFh=findobj('Tag','evAmpFigure','type','figure');
    % initialize figure
    if doMonitorAmpDet
      if isempty(cAmpDetFh)
        cAmpDetFh=figure('Units','normalized', ...
          'Name','PSC Amplitude Plot', ...
          'NumberTitle','off', ...
          'Position',[0.005 0.33 0.99 0.35], ...
          'Tag','evAmpFigure'...
          );
      else
        figure(cAmpDetFh);
        cla
      end
    end
    % lowpass filter as set in threshdetgui (** note that if data had not
    % been lowpass filtered in threshdetgui they won't be filtered here
    % because edge frequencies which are too close to ds.differfi.stopbf
    % are detrimental to the PSC detection algorithm in detPSCAmp.m)
    if ~isempty(head.ap.loCFreq) && isfinite(head.ap.loCFreq)
      if ~isequal(head.ap.loCFreq,loCFreq)
        % overwrite variables
        loCFreq=head.ap.loCFreq;
        loD=lofi(d,si,loCFreq);
      end
    else 
      loD=d;
    end
    % differentiator filter
    if isfield('differfi',head.ap) && ~isempty(head.ap.differfi)
      tmpArgin=eval(head.ap.differfi);
      diffD=differfi(loD,si,tmpArgin{:});  
    else
      % use defaults 
      diffD=differfi(loD,si,ds.differfi.fo,ds.differfi.passbf,ds.differfi.stopbf,...
        'scalFac',ds.differfi.scalFac);
    end
    % convert crucial parameters to points
    tslPts=cont2discrete(evt.tsl{1},si/1e3);
    winEvtCutoutPts=cont2discrete(head.ap.winEvtCutout,si/1e3,'intv',1);
    % create cutouts from data
    evtCutout=tsl2exc(loD,'idx',{tslPts},'win',winEvtCutoutPts);
    % create cutouts from differentated data
    diffEvtCutout=tsl2exc(diffD,'idx',{tslPts},'win',winEvtCutoutPts);
    % call detPSCAmp
    [amp,tRise]=detPSCAmp(evtCutout,diffEvtCutout,1-winEvtCutoutPts(1),...
      si,tslPts,'d',loD,'thresh',head.ap.thresh,'fh',cAmpDetFh,...
      'nPlotEv',min(1000,numel(tslPts))*double(doMonitorAmpDet),...
      'plotOverview',doMonitorAmpDet);
    
    % embed parameters computed above
    for idi=1:nCPscDetPar
      sliceIx=strcmp(ds.cPscDetPar{idi},depPar);
      switch ds.cPscDetPar{idi}
        case 'allCAmp'
          ud.pscr{rowIx,colIx,sliceIx}=amp;
          ud.pscrMn(rowIx,colIx,sliceIx)=nanmedian(amp);
        case 'allCTRise'
          ud.pscr{rowIx,colIx,sliceIx}=tRise(:,1);
          ud.pscrMn(rowIx,colIx,sliceIx)=nanmedian(tRise(:,1));
        case 'allCTRise20_80'
          ud.pscr{rowIx,colIx,sliceIx}=tRise(:,2);
          ud.pscrMn(rowIx,colIx,sliceIx)=nanmedian(tRise(:,2));
        otherwise
          error('bad cPscDet par');
      end
    end
    if doMonitorAmpDet
      switch gr.visualAmpDet
        case 'pause'
          pause(1)
          close(cAmpDetFh)
        case 'waitfor'
          waitfor(cAmpDetFh)
      end
    end
  end
  
  % -----------------------------------------------------------------------
  %          collect/compute phantosic parameters from raw data 
  % -----------------------------------------------------------------------  
  disp('collecting & computing PSC detection-independent parameters from raw data (''phantosic'')...')  
  if ~isempty(ds.loCFreq) && isfinite(ds.loCFreq)
    if ~isequal(ds.loCFreq,loCFreq)
      loCFreq=ds.loCFreq;
      loD=lofi(d,si,loCFreq);
    end
  end
  if any(ismember(depPar,{'baseline','noise','chargePhas'}))
    nPts=ds.phIntv/(si/1e3);
    [base,dev,phas]=phantosic(loD,nPts,ds.phBin,'method',ds.phMethod,...
      'polarity',ds.phPolarity,'prc',ds.phPrc,'frame',gr.doMonitorPhantosic,...
      'pau',1,'si',si);
  end
  % embed phantosic parameters computed above
  for idi=1:nRawPar
    sliceIx=strcmp(ds.rawPar{idi},depPar);
    switch ds.rawPar{idi}
      case 'baseline'
        % baseline for recording
        ud.pscr{rowIx,colIx,sliceIx}=base;
        ud.pscrMn(rowIx,colIx,sliceIx)=nanmean(base);
      case 'noise'
        % deviation from baseline (noise)
        ud.pscr{rowIx,colIx,sliceIx}=dev;
        ud.pscrMn(rowIx,colIx,sliceIx)=nanmean(dev);
      case 'chargePhas'
        % phasic current
        ud.pscr{rowIx,colIx,sliceIx}=phas;
        ud.pscrMn(rowIx,colIx,sliceIx)=nanmean(phas);
      otherwise
        error('bad derived raw par');
    end
  end

  if doBivarHistogram
    disp('generating bivariate histogram...')
    % - transform ud.pscr & extract parameters of interest
    tmppscd=permute(ud.pscr(rowIx,:,:),[1 3 2]);
    pscd=[tmppscd{1,strcmp(ds.plotPar(1).name,depPar),colIx},...
      tmppscd{1,strcmp(ds.plotPar(2).name,depPar),colIx}];
    n=hist3(pscd,'Edges',{ds.plotPar(1).bin2 ds.plotPar(2).bin2});
    % - stuff histogram into appropriate place in ud.histogram
    ud.histogram(:,:,rowIx,colIx)=n;
  end
end
% update userdata and list of experiments
set(masterFh,'userdata',ud);
figure(masterFh);
set(listH,'string',char(ud.listExp));
drawnow

% =========================================================================
%                         PART 3: plot data
% =========================================================================
% currently, if ds.printas is empty, don't plot
if gr.doPlot
  disp('plotting...')
  % restrict to parameters available
  pars=intersect({ds.plotPar.name},depPar,'stable');
  nPars=numel(pars);

  % main figure
  switch gr.ornt
    case {'landscape','portrait'}
      fh1=mkfig(2,'b');
      nPlotCols=3;
      nPlotRows=max(2,ceil(nPars/nPlotCols));
    case 'tall'
      fh1=mkfig(2,'v');
      nPlotCols=2;
      nPlotRows=max(3,ceil(nPars/nPlotCols));
  end
  clf
  orient(gr.ornt);
  
  % 1. -------- univariate histograms of IPSC parameters ---------
  set(gcf,'defaultlinelinewidth',2.5);
  set(gcf,'defaultaxesFontsize',10);
  
  for parIx=1:nPars
    subplot(nPlotRows,nPlotCols,parIx);
    hold on
    % ** rowIx has been computed above and can still be used
    for g=1:nFile
      % index into column of global results variable ud.pscr
      colIx=ds.indepPar==ds.fList{g,2};
      % pull out data
      helpIx=strcmp(depPar,pars{parIx});
      n=ud.pscr{rowIx,colIx,helpIx};
      hh=histc(n,ds.plotPar(parIx).bin);
      plot(ds.plotPar(parIx).bin,hh);
    end
    grid on
    axis tight
    xlabel(ds.plotPar(parIx).name);
    ylabel('N');
  end
  
  legend(ds.indepParLabel,'Interpreter','none');
  
  if ~isempty(gr.printas)
    print(gr.printas,'-r400',[gr.fDir figName '_hist']);
  end
  drawnow
  pause(.05);
  
  
  % 2. ------------ 2D histograms of IPSC parameters -------------
  if doBivarHistogram
    mkfig(3,'b'); clf, orient landscape;
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
    compHist=ud.histogram(:,:,rowIx,compColIx);
    if doNormalize
      % normalize
      compHist=compHist/sum(compHist(:));
    end
    
    sph=gobjects(2,nFile);
    for g=1:nFile
      % index into column of global results variable ud.pscr
      colIx=ds.indepPar==ds.fList{g,2};
      % pull out histogram
      hi=ud.histogram(:,:,rowIx,colIx);
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
      title(ds.indepParLabel{g},'Interpreter','none');
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
    
    if ~isempty(gr.printas)
      print(gr.printas,'-r400',[gr.fDir figName '_2Dhist_' ds.plotPar(1).name '_' ds.plotPar(2).name]);
    end
  end
  
end
drawnow
% focus back on control window
figure(masterFh);


% ------------ LOCAL FUNCTIONS ------------------------------------------
function f=checkFields(f,f_in)
s=fieldnames(f);
s_in=fieldnames(f_in);
sect=setdiff(s_in,s);
if ~isempty(sect)
  errordlg({'nonmatching field names(s): '; char(sect)});
  error('see error dialogue');
end
% assign input values of input struct to current one
for g=1:numel(s_in)
  f.(s_in{g})=f_in.(s_in{g});
end

function loCFreq=determineLoCFreq(ds,varargin)
if nargin>1
  fitHead=varargin{1};
else
  fitHead=[];
end
% default: don't filter
filterCase='';
loCFreq=[];
if ~isempty(fitHead) && isfield(fitHead.ap,'loCFreq') && isfinite(fitHead.ap.loCFreq)
  % data were filtered in pscfitgui
  filterCase=[filterCase '_pscfitgui'];
end
if ~isempty(ds.loCFreq) && isfinite(ds.loCFreq)
  % data were filtered by program in which PSCs were detected
  filterCase=[filterCase '_caller'];
end
switch filterCase
  case '_pscfitgui'
    loCFreq=fitHead.ap.loCFreq;
  case '_caller'
    loCFreq=ds.loCFreq;
  case '_pscfitgui_caller'
    % frequency set in pscfitgui has precedence, but warn
    loCFreq=fitHead.ap.loCFreq;
    if ~isequal(fitHead.ap.loCFreq,ds.loCFreq)
      warndlg({'divergent lowpass corner frequencies:',...
        ['fitHead.ap.loCFreq=' num2str(fitHead.ap.loCFreq) 'Hz'],...
        ['ds.loCFreq=' num2str(ds.loCFreq) 'Hz'],...
        'setting to the former'});
    end
end

function [ds,gr]=pscdeal_defaultParams
  % --- data set descriptors
% the complete set of values of independent parameter (=concentrations)
ds.indepPar=[];
% matching labels
ds.indepParLabel={[]};
% name suffix and extension of files produced by PSCFitgui
ds.fnSuffix='';
% --- analysis parameters
% -- value of independent parameter corresponding to control (against which
% others will be normalized; set to nan for no normalization)
ds.normIpVal=nan;
% -- value of independent parameter of condition against which others will
% be (statistically) COMPARED
ds.compIpVal=[];
% -- parameters computed in pscfitgui (or which can be derived from them)
% to be collected and analyzed
ds.pscFitPar={''};
% -- parameters to be computed from raw data and time stamp lists of PSCs
% using algorithms as used in pscfitgui: the prefix 'all' implies being
% determined from all detected (not necessarily fitted) PSCs, whereas
% parameters lacking this suffix are determined from the subset of detected
% AND fitted PSCs
ds.pscDetPar={''};
% -- parameters to be computed from raw data and time stamp lists of PSCs
% via algorithm geared towards dealing with compound PSCs (occurring in
% bursts with at least partly overlapping rise times):
ds.cPscDetPar={''};
% -- parameters to be determined from raw data only (i.e., independent of
% detected or fitted PSCs)
ds.rawPar={''};
% -- parameters to be plotted (struct)
ds.plotPar=[];
% -- raw data
% corner frequency of lowpass filter (set to [] for no filtering)
ds.loCFreq=nan;
% settings for differentiator filter
% - filter order
ds.differfi.fo=30;
% - lower end of passband (Hz)
ds.differfi.passbf=40;
% - lower end of stopband (Hz, must be substantially higher than
% ds.loCFreq)
ds.differfi.stopbf=3000;
% control PSC amplitude det visualization 
gr.visualAmpDet='none';
% ** all following parameters are input parameters into function
% phantosic.m:
% - set to positive integer if phantosic computations are to be visualized
% (1=each frame, 2=every second frame, ...; 0 for no visuals)
gr.doMonitorPhantosic=0;
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
ds.dDir='';
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
