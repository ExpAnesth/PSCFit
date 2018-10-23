function pscfitguifunc(~,~,job,varargin)
% ** function pscfitguifunc(src,eventdata,job,varargin)
% Collection of callback routines for pscfitgui.m
%                         >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% job              cell array of char    jobs to accomplish
% sp               struct                handles to subplots of main
%                                         figure

% -------------------------------------------------------------------------
% Version 2.4, October 2018
% (C) Harald Hentschke (University Hospital of Tuebingen)
% -------------------------------------------------------------------------

% We need persistent variables:
%   tsl = time stamp list of events
%   tsl_pts = same in discrete time (points)
%   cutout = cutouts of events
%   cutoutFilt = filtered version (for peak detection and plotting purposes)
%   cutoutFitHist = cutouts of 'event history'
%   d = raw data trace
%   dFit = synthetic data trace composed of all PSCS generated from fitted
%   parameters
%   dFitHist = synthetic data trace composed of those PSCS deemed fit for
%   PSC 'history'
% Strucures:
%   ds='data set' listing properties of current file
%   ap='analysis parameters' describing details of current analyis
%   wp='working parameters' (like colors)
%   sp=subplot handles
%   r=results structure
%   head=analysis parameters

persistent ds ap wp sp head r tsl tsl_pts cutout cutoutFilt cutoutFitHist d dFit dFitHist

% to do:
% - variable dFit is not very helpful as it contains all fitted IPSCs
% including the rejected ones

pvpmod(varargin);
etslconst;

done=0;
while ~done
  partJob=job{1};
  switch partJob
    case 'init'
      % *******************************************************************
      % in this ini job, all fields of wp, ap and ds are set up. This is a
      % necessity because the callbacks set here are handles to functions
      % expecting fully defined variables. Besides, it is helpful to have
      % an overview of all fields of above-mentioned key variables. Last,
      % but not least, default values for some parameters are set.
      % *******************************************************************
      % --------------------------
      % ----- set up ds (data set)
      % --------------------------
      % name of major results file
      ds.evtTslFn='lastFile_res';
      % name of file holding evtCutout
      ds.evtCutoutFn='lastFile_evtCutout';
      ds.dataFnCore='lastFile';
      ds.dataPath='e:\_data\_IPSCFit\';
      % +1 for positive-going events, -1 for negative-going events
      ds.polarity=nan;
      % ds.fileInfo=[];

      % ---------------------------------------------------------------
      % ----- set up ap (analysis parameters) with plausible values
      % ---------------------------------------------------------------
      % ~~~~~~~ file names
      % ~~~~~~~ data reading parameters section
      % - min sampling freq (Hz), that is, the data will be downsampled to
      % the nearest value above this which can be achieved by picking
      % interleaved data points
      ap.minSampFreq=1e4;
      % - cutout window (ms)
      ap.winEvtCutout=[-10 40];
      % - pre-trigger dead time (ms) (the posttrigger dead time is
      % equivalent to xIntvFitEndMin below
      ap.preTrigDeadT=0;
      % - lowpass filter freq
      ap.loCFreq=nan;
      % - smoothing of data: span of filter (ms) (span=number of neighbours
      % of a data point that shall make it into the smoothing process)
      ap.smoothTSpan=1;
      % - applies to data smoothing with Savitzky-Golay filters: polynomial
      % order
      ap.smoothPolyOrder=5;
      % ~~~~~~~ data fit parameters section
      % - number of exponential components to fit to current decay
      ap.nDecayComponent=1;
      % - interval in which to compute base line (ms)
      ap.xIntvBaseLine=[-8 -2];
      % - interval in which to seek peak (ms)
      ap.xIntvPeak=[0 4];
      % - minimally acceptable value of right border of interval in which to fit (ms)
      ap.xIntvFitEndMin=10;
      % - maximal value of right border of interval in which to fit (ms)
      ap.xIntvFitEndMax=40;
      % - fitting interval starts postPeakFitDelay ms after peak
      ap.postPeakFitDelay=0.4;
      % - the fraction of maximal PSC amplitude at which width will be
      % assessed (set to 0.5 for half-width)
      ap.widthFac=0.5;
      % upper limit of slow decay component (for biexponential fits) in ms
      ap.maxSlowDecay=100;
      % - was average of cutouts processed?
      ap.isAverageCutout=0;
      % - the minimal fit quality, expressed in terms of adjusted r square,
      % which each event has to fulfil in order to make it into the event
      % history array needed for proper fits to closely spaced PSCs (if set
      % to NaN, history adjustment will not be performed)
      ap.minFitQuality=0.6;
      % - if true, the recent history of the PSCS trace is built into the
      % fitting process, potentially improving proper fits to 'shark fins'
      % (PSCs closely spaced in time)
      ap.doPSCHistory=isfinite(ap.minFitQuality) && ap.minFitQuality<1;
      % - comment on set of parameters
      ap.parComment='a little space for notes';
      
      % §§§ an experimental parameter... tbc
      ap.doGlobalBaseline=false;
      
      % --------------------------------------
      % ----- set up wp ('working' parameters)
      % --------------------------------------
      % ~~~~~~~ display options & matlab version section
      % which version of Matlab?
      wp.mver=ver;
      % note that in standalone deployed code function ver may produce
      % several entries with .Name equal to 'Matlab', so we have to opt for
      % one
      tmpIx=find(strcmpi('matlab',{wp.mver.Name}),1);
      wp.mver=str2double(wp.mver(tmpIx).Version);
      % *** handles to all objects in figure
      wp.handles=guihandles(findobj('tag','PSCFit'));
      % - counterparts to fields of ap but expressed in points
      wp.preTrigDeadT_pts=nan;
      wp.xIntvBaseLine_pts=[nan nan];
      wp.xIntvPeak_pts=[nan nan];
      wp.xIntvFitEndMin_pts=nan;
      wp.xIntvFitEndMax_pts=nan;
      wp.postPeakFitDelay_pts=nan;
      wp.smoothTSpan_pts=nan;
      % interval for detection of 10-90% rise time (will be defined in
      % subfunction convertXInt)
      wp.xIntvRiseT_pts=[nan nan];
      % time axis in pts
      wp.timeAx_pts=nan;
      % offset in pts corresponding to pre-trigger window
      wp.xOffs_pts=nan;
      % y limits of cutouts (mV) - if nan limits will enclose almost all data
      wp.coYLim=[nan nan];
      % max number of cutouts to plot in overlay plot
      wp.maxNPlotCutout=200;
      % standard background color of subplots in main figure window
      wp.stdAxCol=[.85 .85 .85];      
      % color of patch representing x interval for base line
      wp.xIntvBaseLineCol=[.7 .7 .7];      
      % color of patch representing x interval for peak search
      wp.xIntvPeakCol=[1 .5 .1];      
      % color of lines representing right borders of x intervals for fit
      wp.xIntvFitEndMinCol=[.5 .5 1];
      wp.xIntvFitEndMaxCol=[0 0 1];
      % color of fitted events
      wp.FitEventCol=[.9 .8 .2];
      wp.FitEventCol=[.9 .2 .2];
      % the operation mode of cursor in scatter plot
      wp.scatterCursorMode='probe';
      % parameters on display in scatter plot
      wp.scatterD=[];
      % ~~~~~~~ saving options section
      wp.evtTslFnString='';
      % ----- the following working pars are not accessible in the
      % parameters dialog
      % number of rows, columns
      wp.coNRow=0;
      wp.coNCol=0;
      wp.si=nan;
      % x coordinate of last mouse click in scatter plot
      wp.curPtX=0;
      % a flag indicating whether cutouts have been prepared for fit
      wp.isPrepareFit=false;
      % a flag indicating whether cutouts have been fitted
      wp.isFit=false;
      % Hz, hipass corner freq for determination of noise in cutouts
      % (all frequencies above this are supposedly noise)
      wp.noiseHiCFreq=1000;
      % limit of R squared BELOW which events will be kicked out after
      % fitting
      wp.qFitLim=.3;
      % limit of xx ABOVE which events will be kicked out after
      % fitting
      wp.qFit3Lim=1.0;
      % another 'internal' parameter, deciding on whether psc cutouts and
      % related variables are supposed to be dumped on the base workspace
      % after fitting (for illustration and/or monitoring purposes)
      wp.doDumpPSC=true;
      % --------------------------------------
      % ----- set up results struct
      % --------------------------------------
      r=preallocateResults(r,1,1);

      % --------------------------------------
      % ----- initialize subplots
      % --------------------------------------
      % -- cutouts (all)
      % bg and cutout colors 
      set(sp.cutout.axH,'color',wp.stdAxCol,'colororder',coma('bluered','n',50));
      subplot(sp.cutout.axH), hold on,
      % plot the patches representing x intervals
      tmpyl=get(sp.cutout.axH,'ylim');
      sp.cutout.xIntvBaseLinePatchH=patch(wp.xIntvBaseLine_pts([1 1 2 2])',tmpyl([1 2 2 1])',wp.xIntvBaseLineCol);
      sp.cutout.xIntvPeakPatchH=patch(wp.xIntvPeak_pts([1 1 2 2])',tmpyl([1 2 2 1])',wp.xIntvPeakCol);
      % the lines representing minimal and maximal right borders of fit interval
      sp.cutout.xIntvFitEndMinLineH=line([nan nan]',tmpyl([1 2])','color',wp.xIntvFitEndMinCol,'LineWidth',1.5);
      sp.cutout.xIntvFitEndMaxLineH=line([nan nan]',tmpyl([1 2])','color',wp.xIntvFitEndMaxCol,'LineWidth',1.5);
      % handle to cutouts
      sp.cutout.ph=plot(nan);
      % handle to text
      sp.cutout.th=smarttext('no data loaded',...
        .99,.03,'color','k','fontsize',10,'fontweight','bold','interpreter','none');
      % index to cutouts currently marked in scatter plot
      wp.curScatCoIx=[];
      % handle to those cutous in scatter plot
      wp.curScatCoPh=[];
      % bg color for selected cutouts plot
      set(sp.selectCutout.axH,'color',wp.stdAxCol);
      % bg color for scatter plot
      set(sp.scatter.axH,'color',wp.stdAxCol);      
      % next two jobs: set scatter cursor mode and set values of uicontrols
      job(2:end+1)=job;
      job(1:2)={'setScatterCursorMode','writeParameters2Gui'};

    case 'loadParametersFromFile'
      % see job 'writeParameters2File' for comprehensive explanation
      OKFlag=0;
      uicFn=fieldnames(wp.handles);
      % apFn contains names of all fields that can be set by uicontrols
      apFn=intersect(fieldnames(ap),uicFn);
      % same for wp
      wpFn=intersect(fieldnames(wp),uicFn);
      if ~isdeployed
        % by default look for files in \PSCFit\parameterFiles
        pfDir=mfilename('fullpath');
        pfDir=pfDir(1:max(strfind(pfDir,'\'))-1);
        w=what;
        cd([pfDir '\parameterFiles']);
      end
      [tmpOptFn,tmpOptPath] = uigetfile('*.mat','pick parameter file');
      if ischar(tmpOptFn) && ischar(tmpOptPath)
        load([tmpOptPath tmpOptFn],'ap_uic','wp_uic');
        OKFlag=1;
        if ~exist('ap_uic','var') || ~exist('wp_uic','var')
          errordlg('Chosen file does not contain parameters. Please choose a different file.');
          OKFlag=0;
        end
        % *** make sure that fields of ap_uic and wp_uic match exactly with
        % the ui-controllable fields of ap and wp (as defined above
        % (apFn=...))
        if ~isempty(setxor(fieldnames(ap_uic),apFn)) || ~isempty(setxor(fieldnames(wp_uic),wpFn))
          errordlg('illegal or missing parameter(s) in file - possibly an outdated parameter file was loaded. Please choose a different file.');
          OKFlag=0;
        end
        if OKFlag
          % assign values to fields of ap and wp
          for g=1:numel(apFn)
            ap.(apFn{g})=ap_uic.(apFn{g});
          end
          for g=1:numel(wpFn)
            wp.(wpFn{g})=wp_uic.(wpFn{g});
          end
        end
      end
      if ~isdeployed
        % cd back to original dir
        cd(w.path);
      end
      if OKFlag
        % next two jobs: check parameters and write them to gui
        job(2:end+1)=job;
        job(1:2)={'digestParameters','writeParameters2Gui'};
      else
        job(1)=[];
      end
      clear w pfDir ap_uic wp_uic OKFlag

    case 'writeParameters2Gui'
      % set the various uicontrols' strings and values to those of
      % corresponding fields of ap and wp, all along checking for errors
      uicFn=fieldnames(wp.handles);
      apFn=fieldnames(ap);
      wpFn=fieldnames(wp);
      structName={'ap','wp'};
      % set 'string' properties of handles to all uicontrols to the values
      % of the matching fields of ap or wp
      for g=1:length(uicFn)
        structIx=[any(strcmp(uicFn{g},apFn)),...
          any(strcmp(uicFn{g},wpFn))];
        if length(find(structIx))==1
          eval(['cType=get(wp.handles.' uicFn{g} ',''style'');']);
          switch cType
            case 'edit'
              % depending on the type of the field...
              switch uicFn{g}
                case {'parComment','evtTslFnString'}
                  eval(['set(wp.handles.' uicFn{g} ',''string'',' structName{structIx} '.' uicFn{g} ');']);
                  %                 case {'thresh'}
                  %                   % threshold must be written back with comparatively high
                  %                   % precision
                  %                   eval(['set(wp.handles.' uicFn{g} ',''string'',num2str(' structName{structIx}  '.' uicFn{g} ',''% 2.5f''));']);
                  %                 case {'excYLim'}
                  %                   % same for y limits of excerpt plot
                  %                   eval(['set(wp.handles.' uicFn{g} ',''string'',num2str(' structName{structIx}  '.' uicFn{g} ',''% 4.3f''));']);
                otherwise
                  eval(['set(wp.handles.' uicFn{g} ',''string'',num2str(' structName{structIx}  '.' uicFn{g} ',''% 8.2f''));']);
              end
            case 'checkbox'
              eval(['set(wp.handles.' uicFn{g} ',''value'',' structName{structIx}  '.' uicFn{g} ');']);
            case 'text'
              % do nothing because text uicontrols do not hold any
              % information
            otherwise
              errordlg('internal: encountered uicontrol other than edit, checkbox or text');
          end
        elseif length(find(structIx))>1
          errordlg(['internal: uicontrol tagged ''' uicFn{g} ''' has more than one matching fields in ap and wp']);
        end
      end
      job(1)=[];

    case 'readParametersFromGui'
      % the inverse of job 'writeParameters2Gui': retrieve the various
      % uicontrols' strings and values and set corresponding fields of ap
      % and wp. All checks for major pitfalls are done in 'writeParameters2Gui'
      % so they are omitted here.
      uicFn=fieldnames(wp.handles);
      apFn=fieldnames(ap);
      wpFn=fieldnames(wp);
      structName={'ap','wp'};
      for g=1:length(uicFn)
        structIx=[any(strcmp(uicFn{g},apFn)),...
          any(strcmp(uicFn{g},wpFn))];
        if length(find(structIx))==1
          eval(['cType=get(wp.handles.' uicFn{g} ',''style'');']);
          switch cType
            case 'edit'
              % depending on the type of the field...
              switch uicFn{g}
                case {'parComment','evtTslFnString'}
                  eval([structName{structIx} '.' uicFn{g} '=get(wp.handles.' uicFn{g} ',''string'');']);
                otherwise
                  eval(['[tmpNum,conversionOK]=str2num(get(wp.handles.' uicFn{g} ',''string''));']);
                  if conversionOK
                    eval([structName{structIx} '.' uicFn{g} '=tmpNum;']);
                  else
                    warndlg(['could not read value for ' structName{structIx} '.' uicFn{g} ' - only numeric values are allowed (typographic error?)'])
                  end
              end
            case 'checkbox'
              eval([structName{structIx}  '.' uicFn{g} '=get(wp.handles.' uicFn{g} ',''value'');']);
            otherwise
          end
        end
      end
      job(1)={'digestParameters'};

    case 'saveParameters2File'
      % *** there is no 1:1 correspondence between the uicontrols and
      % fields of ap and wp. For example, evtTslFn is a field of ap that
      % cannot be set via a uicontrol but depends on the specific file that
      % is currently analyzed. Such file-specific information must not be
      % saved. Instead, only those fields of ap and wp represented by
      % uicontrols shall be saved. This is done in the lines below
      uicFn=fieldnames(wp.handles);
      % ap_uic contains all fields that can be set by uicontrols
      ap_uic=rmfield(ap,setdiff(fieldnames(ap),uicFn));
      % same for wp
      wp_uic=rmfield(wp,setdiff(fieldnames(wp),uicFn));
      if ~isdeployed
        % by default dump files in \PSCFit\parameterFiles
        pfDir=mfilename('fullpath');
        pfDir=pfDir(1:max(strfind(pfDir,'\'))-1);
        w=what;
        cd([pfDir '\parameterFiles']);
      end
      [tmpDataFn,tmpDataPath] = uiputfile('*.mat');
      if ischar(tmpDataFn) && ischar(tmpDataPath)
        save([tmpDataPath tmpDataFn],'ap_uic','wp_uic');
      end
      if ~isdeployed
        % cd back to original dir
        cd(w.path);
      end
      % delete vars
      clear w pfDir ap_uic wp_uic
      job(1)=[];

    case 'digestParameters'
      % this part job is run whenever parameters were read from gui
      job(1)=[];
      isAllParameterOK=true;
      doDigest=true;
      % previous results must be discarded
      [wp,r]=discardResults(wp,sp);
      if doDigest
        disp('** processing & checking options..');
        % ----- checks of parameters:
        % --- sampling freq
        if ap.minSampFreq<0 || ~isfinite(ap.minSampFreq)
          warndlg('minimal sampling freq must be finite and >0 - setting to 10000');
          isAllParameterOK=false;
          ap.minSampFreq=1e4;
        end
        % --- lopass filter
        if isfinite(ap.loCFreq) && ap.loCFreq<0
          warndlg('lowpass filter cutoff frequency must be a positive value or NaN (=no filter)');
          isAllParameterOK=false;
          ap.loCFreq=nan;
        end        
        % --- smoothing parameter
        if ~isfinite(ap.smoothTSpan)
          warndlg('smoothing time span must be a finite, positive value');
          isAllParameterOK=false;
        end          
        if ap.smoothTSpan<=ap.smoothPolyOrder*1000/ap.minSampFreq
          warndlg('smoothing time span is too small for filtering - if data shall be filtered, increase the value');
          % ** user has been warned, but do not force filtering
        end
        % --- number of decay components
        if isempty(intersect(ap.nDecayComponent, [1 2]))
          warndlg('number of decay components must be 1 or 2 - setting to 1');
          isAllParameterOK=false;
          ap.nDecayComponent=1;
        end
        % --- x intervals
        % - cutout interval
        if length(ap.winEvtCutout)~=2 || any(~isfinite(ap.winEvtCutout))
          warndlg('cutout interval must contain two finite, non-nan values')
          isAllParameterOK=false;
          ap.winEvtCutout=[nan nan];
        elseif diff(ap.winEvtCutout)<=0
          warndlg('left value of base line interval must be smaller than the right one');
          isAllParameterOK=false;
          ap.winEvtCutout=[nan nan];
        end
        % - dead time
        if length(ap.preTrigDeadT)~=1 || ~isfinite(ap.preTrigDeadT)
          warndlg('pre-trigger dead time must be one finite, non-nan value')
          isAllParameterOK=false;
          ap.preTrigDeadT=nan;
        end
        % - base line
        if length(ap.xIntvBaseLine)~=2 || any(~isfinite(ap.xIntvBaseLine))
          warndlg('base line interval must contain two finite, non-nan values')
          isAllParameterOK=false;
          ap.xIntvBaseLine=[nan nan];
        elseif diff(ap.xIntvBaseLine)<=0
          warndlg('left value of base line interval must be smaller than the right one');
          isAllParameterOK=false;
          ap.xIntvBaseLine=[nan nan];
        end
        % - peak search
        if length(ap.xIntvPeak)~=2 || any(~isfinite(ap.xIntvPeak))
          warndlg('peak search interval must contain two finite, non-nan values')
          isAllParameterOK=false;
          ap.xIntvPeak=[nan nan];
        elseif diff(ap.xIntvPeak)<=0
          warndlg('left value of peak search interval must be smaller than the right one');
          isAllParameterOK=false;
          ap.xIntvPeak=[nan nan];
        end
        % - right borders of minimal/maximal fit intervals
        if length(ap.xIntvFitEndMin)~=1 || ~isfinite(ap.xIntvFitEndMin)
          warndlg('right border of minimal fitting interval must contain one finite, non-nan value')
          isAllParameterOK=false;
          ap.xIntvFitEndMin=nan;
        end
        if length(ap.xIntvFitEndMax)~=1 || ~isfinite(ap.xIntvFitEndMax)
          warndlg('right border of maximal fitting interval must contain one finite, non-nan value')
          isAllParameterOK=false;
          ap.xIntvFitEndMax=nan;
        end
        if ap.preTrigDeadT>0
          warndlg('pre-trigger dead time must be zero or negative');
          isAllParameterOK=false;
          ap.preTrigDeadT=nan;
        end
        if diff([ap.xIntvFitEndMin  ap.xIntvFitEndMax])<0
          warndlg('right border of maximal fitting interval must be equal to or more positive than that of minimal fitting interval');
          isAllParameterOK=false;
          ap.xIntvFitEndMax=nan;
        end
        if diff([ap.xIntvPeak(2) ap.xIntvFitEndMin])<=0
          warndlg('right border of minimal fitting interval must be outside of (=more positive than) the peak search interval');
          isAllParameterOK=false;
          ap.xIntvFitEndMin=nan;
        end
        % - post-peak fit delay
        if length(ap.postPeakFitDelay)~=1 || ~isfinite(ap.postPeakFitDelay) || ap.postPeakFitDelay<0
          warndlg('post-peak fit delay must contain one positive, finite, non-nan value')
          isAllParameterOK=false;
          ap.postPeakFitDelay=nan;
        end
        if diff([ap.postPeakFitDelay  ap.xIntvFitEndMin])<0
          warndlg('post-peak fit delay is far too large');
          isAllParameterOK=false;
          ap.postPeakFitDelay=nan;
        end
        % - required minimal fit quality for PSC history-based adjustment
        % of PSC waveform
        if isfinite(ap.minFitQuality)
          ap.doPSCHistory=true;
          if ap.minFitQuality>1
            warndlg('minimal fit quality for PSC waveform adjustment must be smaller than 1.0 - setting to NaN');
            isAllParameterOK=false;
            ap.minFitQuality=NaN;
          elseif ap.minFitQuality>.95
            warndlg('minimal fit quality for PSC waveform adjustment is very high - currently chosen value will be kept, but you''d better make sure this is what you intend');
          elseif ap.minFitQuality<.2
            warndlg('minimal fit quality for PSC waveform adjustment is very low - currently chosen value will be kept, but you''d better make sure this is what you intend');
          end
        else
          ap.doPSCHistory=false;
        end
      end
      % if data is loaded
      % - adjust the values expressed in pts
      % - check whether values are compatible with cutout length
      % - if so, adjust patches
      if isfinite(wp.si)
        wp=convertXInt(ds,ap,wp);
        if wp.xIntvFitEndMax_pts > wp.timeAx_pts(end)
          warndlg('right border of maximal fitting interval is beyond cutout length')
          isAllParameterOK=false;
          ap.xIntvFitEndMax=nan;
        end
        % - no need to check xIntvFitEndMin, postPeakFitDelay or
        % preTrigDeadT_pts
        if wp.xIntvBaseLine_pts(end) > wp.timeAx_pts(end)
          warndlg('right border of base line interval is beyond cutout length')
          isAllParameterOK=false;
          ap.xIntvBaseLine=[nan nan];
        end
        if wp.xIntvBaseLine_pts(1) < wp.timeAx_pts(1)
          warndlg('left border of base line interval is beyond cutout length')
          isAllParameterOK=false;
          ap.xIntvBaseLine=[nan nan];
        end
        if wp.xIntvPeak_pts(end) > wp.timeAx_pts(end)
          warndlg('right border of peak search interval is beyond cutout length')
          isAllParameterOK=false;
          ap.xIntvPeak=[nan nan];
        end
        if wp.xIntvPeak_pts(1) < wp.timeAx_pts(1)
          warndlg('left border of peak search interval is beyond cutout length')
          isAllParameterOK=false;
          ap.xIntvPeak=[nan nan];
        end
        if isAllParameterOK
          adjustXIntPatch(sp,wp);
        end
      end

      if ~isAllParameterOK
        job(1)={'writeParameters2Gui'};
      end

    case {'clearData','readData'}
      % ******************************************************************
      % re-enable minSampFreq and winEvtCutout uicontrols
      set(wp.handles.minSampFreq,'enable','on');
      set(wp.handles.winEvtCutout,'enable','on');
      set(wp.handles.loCFreq,'enable','on');
      % delete/reset results/plots/variables from last file/channel
      d=[];
      dFit=[];
      dFitHist=[];
      cutout=[];
      cutoutFilt=[];
      cutoutFitHist=[];
      wp.si=nan;
      tsl=[];
      wp.timeAx_pts=nan;
      ds.polarity=nan;
      ap.isAverageCutout=0;
      % this takes care of all fit results, wp.is... and plots, except
      % cutouts plot...
      [wp,r]=discardResults(wp,sp);
      % ... done here
      if any(ishandle(sp.cutout.ph))
        delete(sp.cutout.ph);
        sp.cutout.ph=nan;
      end
      set(sp.cutout.th,'string','no data loaded');
      if strcmp(partJob,'readData')
        doLoad=true;
        [tmpDataFn,tmpDataPath] = uigetfile('*.mat','pick *res* file',ds.dataPath);
        if isnumeric(tmpDataFn) && ~tmpDataFn
          doLoad=false;
        else
          % ** note that only *res* files will be accepted
          if ischar(tmpDataFn) && ~contains(tmpDataFn,'res')
            warndlg('please pick a *res*.mat file');
            doLoad=false;
          end
        end
        if doLoad
          % names of .mat files:
          tmpix=strfind(tmpDataFn,'_res.mat');
          % - 'core' file name
          ds.dataFnCore=tmpDataFn(1:tmpix(end)-1);
          % - path
          ds.dataPath=tmpDataPath;
          % - full name of results file picked above
          ds.evtTslFn=[tmpDataPath ds.dataFnCore '_res'];
          % - full name of accompanying cutouts file (if any)
          ds.evtCutoutFn=[tmpDataPath ds.dataFnCore '_evtCutout'];
          % load tsl and head and detach
          load(ds.evtTslFn,'evt','head');
          % evt.tsl may be an empty array, so catch this
          if iscell(evt.tsl)
            tsl=evt.tsl{1};
          else
            tsl=[];
          end
          % load raw data
          [d,si]=abfload([tmpDataPath head.ds.dataFn],'channels',head.wp.dataChanName);
          % filter
          if isfinite(ap.loCFreq)
            d=lofi(d,si,ap.loCFreq);
          end
          % if downsampling is requested and possible by taking interleaved
          % points, do it now
          tmpSi=1e6/ap.minSampFreq;
          tmpSampFac=floor(tmpSi/si);
          if tmpSampFac>1
            d=d(1:tmpSampFac:end,:);
            si=si*tmpSampFac;
            disp(['data were downsampled to ' num2str(1e6/si) ' Hz']);
          end
          % copy sampling interval to wp
          wp.si=si;
          % continuous traces in which fits of events will be placed
          dFit=zeros(size(d));
          dFitHist=zeros(size(d));
          % length of d
          wp.nd=numel(d);
          % disable minSampFreq and winEvtCutout uicontrols
          set(wp.handles.minSampFreq,'enable','off');
          set(wp.handles.winEvtCutout,'enable','off');
          set(wp.handles.loCFreq,'enable','off');
          % next job: generate cutouts
          job(1)={'genCutout'};
          clear tmp* bu evt
        else
          job(1)=[];
        end
      else
        job(1)=[];
      end
      
    case 'genCutout'
      if ap.doGlobalBaseline
        % §§§§§
        d=detrend(d);
        % §§ should do a much better job here, possibly piecewise linear
        % detrending procedure, picking stretches of data with few events
      end
      % generate cutouts
      [cutout,isCutout]=tsl2exc(d,wp.si,{tsl},'win',ap.winEvtCutout);
      % 'no event left behind (?)': as the current cutout window may differ
      % from the original one make sure that lost cutouts are also removed
      % from etsl
      tsl=tsl(isCutout{1});
      % size of things
      [wp.coNRow,wp.coNCol]=size(cutout);
      % now convert time intervals to pts
      wp=convertXInt(ds,ap,wp);
      % convert time stamps to pts
      tsl_pts=cont2discrete(tsl,wp.si/1000,'intv',0);
      % polarity: if 
      % [average of points in peak search interval - average of points in baseline interval]
      % is positive, we're dealing with positive-going events, and vice versa
      ixPeak=(wp.xIntvPeak_pts(1):wp.xIntvPeak_pts(2))-wp.xOffs_pts;
      ixBase=(wp.xIntvBaseLine_pts(1):wp.xIntvBaseLine_pts(2))-wp.xOffs_pts;
      ds.polarity=sign(mean(median(cutout(ixPeak,:))-median(cutout(ixBase,:))));
      % **************************************************************
      % invert negative-going events to streamline all following
      % operations
      % **************************************************************
      if ds.polarity<0
        cutout=cutout*-1;
      end
      % next job: plot data
      job(1)={'plotOv'};

    case 'plotOv'
      % if this partjob is requested, clear excerpts, if any
      if any(ishandle(sp.cutout.ph))
        delete(sp.cutout.ph);
        sp.cutout.ph=nan;
      end
      if ~isempty(cutout)
        % plot specified number of cutouts
        tmpv1=min(wp.maxNPlotCutout,wp.coNCol);
        tmpix=unique(ceil((1:tmpv1)/tmpv1*wp.coNCol));
        subplot(sp.cutout.axH)
        % - set x and y limits of axis to new values
        set(sp.cutout.axH,...
          'xlim',wp.timeAx_pts([1 end]),...
          'ylim',[min(min(cutout(:,tmpix)))  max(max(cutout(:,tmpix)))],'colororder',coma('bluered','n',50));
        sp.cutout.ph=plot(wp.timeAx_pts,cutout(:,tmpix));
        % set real time units as labels (such that zero and integer values
        % are shown. This can be achieved easily by shifting x ticks as
        % defined by wp.timeAx_pts one tick to the right)
        tmpxtick=cont2discrete(5,wp.si/1000,'intv',1);
        tmpxtick=[fliplr(1:-tmpxtick:wp.timeAx_pts(1)) 1+tmpxtick:tmpxtick:wp.timeAx_pts(end)];
        set(sp.cutout.axH,'xtick',tmpxtick);
        set(sp.cutout.axH,'xticklabel',discrete2cont(tmpxtick,wp.si/1000,'intv',0));
        adjustXIntPatch(sp,wp);
        % information on file
        if ishandle(sp.cutout.th)
          delete(sp.cutout.th);
        end
        sp.cutout.th=smarttext([ds.dataFnCore ', si=' int2str(wp.si) ' us, ' int2str(tmpv1) ' of ' int2str(wp.coNCol) ' events'],...
          .99,.03,'color',[1 1 .4],'fontsize',14,'fontweight','bold','interpreter','none');
        if ap.isAverageCutout
          set(sp.cutout.th,'string',[ds.dataFnCore ', AVERAGE']);
        end
      end
      job(1)=[];
      
    case 'toggleCutoutLineStyle'
      if any(ishandle(sp.cutout.ph))
        if strcmpi(get(sp.cutout.ph,'marker'),'none')
          set(sp.cutout.ph,'marker','.');
        else
          set(sp.cutout.ph,'marker','none');
        end
      end
      job(1)=[];
              
    case 'averageCutout'
      if ~isempty(cutout)
        if ap.isAverageCutout
          warndlg('cutouts have already been averaged');
          uiwait(h);
          job(1)=[];
        else
          % delete/reset results/plots/variables from last file/channel
          cutoutFilt=[];
          cutoutFitHist=[];
          % this takes care of all fit results and plots except sp.cutout.axH
          [wp,r]=discardResults(wp,sp);
          ap.isAverageCutout=1;
          % *** now compute average of cutouts and set ts to 0
          cutout=mean(cutout,2);
          tsl=0;
          tsl_pts=0;
          % size of things
          [wp.coNRow,wp.coNCol]=size(cutout);        
          % next job: re-plot data
          job(1)={'plotOv'};
        end
      else
        h=warndlg('no data loaded');
        uiwait(h);
        job(1)=[];
      end

    case 'prepareFit'
      if ~isempty(cutout)
        % discard all previous results
        [wp,r]=discardResults(wp,sp);
        % set flag that data have been prepared
        wp.isPrepareFit=true;
        % reset continuous traces of fit (history)
        dFit=zeros(size(d));
        dFitHist=zeros(size(d));
        % preallocation
        r=preallocateResults(r,wp.coNCol,ap.nDecayComponent);
        
        % i) base line
        if ap.doGlobalBaseline
          % do NOT subtract base line
        else
          ix=(wp.xIntvBaseLine_pts(1):wp.xIntvBaseLine_pts(2))-wp.xOffs_pts;
          % subtract from traces, don't keep value
          cutout=cutout-repmat(mean(cutout(ix,:)),wp.coNRow,1);
        end
        
        % ii) peaks, following Savitzky-Golay filtering (§§ could be a
        % little more flexible by allowing poly order to go down to 3 or
        % even 2)
        if wp.smoothTSpan_pts>ap.smoothPolyOrder
          cutoutFilt=sgolayfilt(cutout,ap.smoothPolyOrder,wp.smoothTSpan_pts);
        else
          warning(['not filtering traces because smoothing time span is too short (must be >= ' int2str(ap.smoothPolyOrder) ' points']);
          cutoutFilt=cutout;
        end
        % ** note that tPeak is in pts, corrected for xIntvPeak interval
        % such that continuous time zero corresponds to discrete time 1, as
        % all _pts parameters
        ix=(wp.xIntvPeak_pts(1):wp.xIntvPeak_pts(2))-wp.xOffs_pts;
        % for each PSC individually, set ignore all points which occur
        % starting 'a little before' the next detected PSC (if any):
        
        % - this is 'a little before', which is also the minimal number of
        % points that must remain of the cutout for peak detection
        buffer_pts=cont2discrete(0.5,wp.si/1000);
        % - local copy of cutoutFilt
        tmpD=cutoutFilt(ix,:);
        % - ISI
        tmpISI_pts=diff(tsl_pts);
        for ii=1:size(cutoutFilt,2)-1
          if tmpISI_pts(ii)-buffer_pts<=wp.xIntvPeak_pts(2)
            % set to zero so that any subsequent peak will be out of the
            % game
            tmpD(max(tmpISI_pts(ii)-wp.xIntvPeak_pts(1)-buffer_pts,buffer_pts):end,ii)=0;
          end
        end
        evr=evdeal(tmpD,'idx',{'minmaxpeak'});
        r.aPeak=evr.maxPeak;
        r.tPeak=evr.maxPeakT+ix(1)-1+wp.xOffs_pts;
        clear tmpD tmpISI
        
        % iii) 10-90% rise time (determine from SG-filtered version, too)
        ix=(wp.xIntvRiseT_pts(1):wp.xIntvRiseT_pts(2))-wp.xOffs_pts;
        % upsample waveforms to 100 kHz (si=10 us) for more precision (edge
        % effects of built-in filter of function resample should not play a
        % role as baseline is subtracted already and the right edge, prone
        % to distortion, is far beyond the peak)
        iSi=10;
        iCutout=resample(cutoutFilt(ix,:),wp.si,iSi);
        niSample=size(iCutout,1);
        [~,i1]=max(iCutout>=repmat(.1*r.aPeak,niSample,1));
        [~,i2]=max(iCutout>=repmat(.9*r.aPeak,niSample,1));
        r.tRise=(i2-i1)*(iSi/wp.si);
        % assign values to fields of r, adjusting time frame as above
        r.t10=i1*(iSi/wp.si)+ix(1)-1+wp.xOffs_pts;
        r.t90=i2*(iSi/wp.si)+ix(1)-1+wp.xOffs_pts;
        
        % iv) tag all events which are problematic:
        % - no peak found
        r.OKIx=isfinite(r.aPeak);
        disp([int2str(numel(find(r.OKIx))) ' of ' int2str(wp.coNCol) ' events with proper peak']);
        % - has a neighbor somewhere in the interval
        % [wp.preTrigDeadT_pts       wp.xIntvFitEndMin_pts]
        % -- deal with preceding events (watch the sign of things!)
        tmpOKIx=[true; diff(tsl_pts)> -wp.preTrigDeadT_pts];
        % -- deal with lagging events
        tmpOKIx=tmpOKIx & [diff(tsl_pts)>wp.xIntvFitEndMin_pts; true];
        r.OKIx=r.OKIx & tmpOKIx';
        disp([int2str(numel(find(r.OKIx))) ' of ' int2str(wp.coNCol) ' events with proper peak and no immediate neighbors']);

        % v) plot a selection of tmpN cutouts in raw and filtered versions
        % so we can judge whether filters pars are appropriate. Take events
        % from whole pool of cutouts, not just the ones deemed OK, those
        % representing the range of amplitudes but avoiding extremes
        tmpN=16;
        [~,tmpCoIx]=sort(abs(r.aPeak));
        tmpCoIx=tmpCoIx(round(linspace(1,wp.coNCol,tmpN+2)));
        tmpCoIx=unique(tmpCoIx(2:end-1));
        plotProcessTrace(wp,r,cutout,cutoutFilt,sp,tmpCoIx);
      else
        h=warndlg('no data loaded');
        uiwait(h);
        wp.isPrepareFit=false;
      end
      job(1)=[];

    case 'fit'
      % check whether job 'prepareFit' has run before
      if wp.isPrepareFit
        wp.isFit=true;
        % reset continuous traces of fit (history)
        dFit=zeros(size(d));
        dFitHist=zeros(size(d));
        % wipe plots
        cla(sp.selectCutout.axH);
        cla([sp.scatter.axH])
        wp.scatterD=[];
        wp.curScatCoPh=[];
        wp.curScatCoIx=[];
        switch ap.nDecayComponent
          case 1
            ft=fittype('a1*exp(-t/tau1)',...
              'independent',{'t'},...
              'coefficients',{'a1','tau1'});
            % upper limits: make plausibility assumptions
            % - amplitude cannot be more than x*maximal peak found
            % - decay time is extremely unlikely to be more than
            % 2*posttrigger cutout length
            u=[max(r.aPeak)*2  wp.timeAx_pts(end)*2];
            % lower bounds: 
            l=[0 0];
            % anonymous function to flesh out fit with real numbers
            % (currently not needed)
            % fifu=@(t,a1,tau1) a1*exp(-t/tau1);
          case 2
            ft=fittype('a1*exp(-t/tau1) + a2*exp(-t/tau2)',...
              'independent',{'t'},...
              'coefficients',{'a1','tau1','a2','tau2'});
            % upper limits: make plausibility assumptions
            % - either amplitude cannot be more than x*maximal peak found
            % - fast decay time is extremely unlikely to be more than
            % 2*posttrigger cutout length, slow decay time is fixed at
            % wp.maxSlowDecayPts
            u=[max(r.aPeak)*2  wp.timeAx_pts(end)*2  max(r.aPeak)*2  wp.maxSlowDecay_pts];
            % lower bounds: 
            l=[0 0 0 0];
            % anonymous function to flesh out fit with real numbers
            % (currently not needed)
            % fifu=@(t,a1,tau1,a2,tau2) a1*exp(-t/tau1) + a2*exp(-t/tau2);
        end
        
        fo=fitoptions('method','NonlinearLeastSquares',...
          'upper',u,...
          'lower',l,...
          'TolX',.001,...
          'display', 'off'...
          );

          % 'robust','Bisquare',...
        
        % define right border of fit interval as described in detail below
        % and store this information in field of r (keep distance of
        % wp.xIntvPeak_pts/4 to next event in order to prevent its rise
        % phase influencing the present one)
        r.xIntvFitEnd=min(wp.xIntvFitEndMax_pts,[(diff(tsl_pts)' -round(diff(wp.xIntvPeak_pts)/4)) inf]);
        
        % generate hipass-filtered version of cutouts from which base line
        % noise magnitude can be computed (will be needed further down to
        % produce one of three measures of the goodness of the fit)
        hifiCutout=hifi(cutout,wp.si,wp.noiseHiCFreq).^2;
        
        % container for fits to cutouts
        cutoutFit=zeros(size(cutout));
        
        % also preallocate container for event history
        cutoutFitHist=zeros(size(cutout));
        
        % recompute index to base lines of events
        baseIx=(wp.xIntvBaseLine_pts(1):wp.xIntvBaseLine_pts(2))-wp.xOffs_pts;
        % set base line to zero for all events
        r.base=zeros(size(r.base));
        
        h=waitbar(0,'Fitting...');
        for g=1:wp.coNCol
          if r.OKIx(g)
            % first, the interval to be fitted must be determined:
            % - the first point to be fitted is wp.postPeakFitDelay_pts
            % after the peak
            % - the last point is, ideally, xIntvFitEndMax_pts
            % - ix3 is the time axis (unit: points) extending to this
            % point, used for generating fits to the cutouts and
            % determination of the (half-) width
            % - if however, there is a following event in this interval,
            % move the interval border back towards zero, but only up to
            % xIntvFitEndMin_pts (all events with a closer event have been
            % kicked out in job preparefit): this is ix2
            % time axis for fit: unit is still ticks, but we start at zero
            ix=(0:r.xIntvFitEnd(g)-r.tPeak(g)-wp.postPeakFitDelay_pts)';
            % the corresponding index into cutout
            ix2=(r.tPeak(g)+wp.postPeakFitDelay_pts:r.xIntvFitEnd(g))'-wp.xOffs_pts;
            % time axis for generating fits based on parameters extracted
            % from fit
            ix3=(0:wp.xIntvFitEndMax_pts-r.tPeak(g)-wp.postPeakFitDelay_pts)';
            % the corresponding index into cutout
            ix4=(r.tPeak(g)+wp.postPeakFitDelay_pts:wp.xIntvFitEndMax_pts)'-wp.xOffs_pts;
            if ap.doPSCHistory && g>1 && r.qFit(g-1)>ap.minFitQuality
              % 1. cut out history spanning entire present cutout and place
              % in cutoutFitHist
              curHistCutout=dFitHist(tsl_pts(g)+wp.timeAx_pts);
              cutoutFitHist(:,g)=curHistCutout;
              % 2. subtract it from present full cutout
              curCutout=cutout(:,g)-curHistCutout;
              % 3. re-adjust base line, keep the value & cut out decay phase
              r.base(g)=mean(curCutout(baseIx));
              curDecayCutout=curCutout(ix2)-r.base(g);
            else              
              % decay phase of present cutout
              curDecayCutout=cutout(ix2,g);
            end
            
            try
              if ap.nDecayComponent==1
                % start values for fit: deduce from parameters computed
                % thus far
                set(fo,'Startpoint',[r.aPeak(g) wp.xIntvFitEndMin_pts]);
                [f,gof]=fit(ix,curDecayCutout,ft,fo);
                % extract parameters
                tmpPar=coeffvalues(f);
                % amplitude of decay component
                r.amp(1,g)=tmpPar(1);
                % decay time(s)
                r.tDecay(1,g)=tmpPar(2);
                % determine width of fitted function:
                % - fit
                curCutoutFit=f(ix3);
                % - indexes to points of fit flanking intersection with zero 
                x1=find(curCutoutFit>r.amp(end,g)*ap.widthFac,1,'last');
                if ~isempty(x1) && x1<=ix(end)
                  x2=x1+1;
                  % corresponding amplitudes
                  y1=curCutoutFit(x1);
                  y2=curCutoutFit(x2);
                  % intersection:
                  r.width(1,g)=x1+(r.amp(end,g)*ap.widthFac-y1)*(x2-x1)/(y2-y1);
                end
                % now embed
                cutoutFit(ix4,g)=curCutoutFit;
                % quality of fit:
                % 1. adjusted r square
                r.qFit(g)=gof.adjrsquare;
                % 2. root mean square error
                r.qFit2(g)=gof.rmse;
                % 3. mean squared error minus mean squared error of noise
                % (as obtained via hipass filtering of traces), normalized
                % to PSC's amplitude
                r.qFit3(g)=(gof.sse - sum(hifiCutout(ix2,g)))/numel(ix)/r.amp(1,g);
              else
                % by default, assume that fast component has higher amplitude
                set(fo,'Startpoint',[r.aPeak(g)*.7  wp.xIntvFitEndMin_pts/2   r.aPeak(g)*.3  wp.xIntvFitEndMin_pts*2]);
                [f,gof]=fit(ix,curDecayCutout,ft,fo);
                tmpPar=coeffvalues(f);
                % in case the fit converged to tau1 as the slower component
                % switch values
                if diff(tmpPar([2 4]))<0
                  tmpPar=tmpPar([3 4 1 2]);
                end
                % amplitudes of decay component(s)
                r.amp(1,g)=tmpPar(1);
                r.amp(2,g)=tmpPar(3);
                % decay time(s)
                r.tDecay(1,g)=tmpPar(2);
                r.tDecay(2,g)=tmpPar(4);
                % summed amplitudes
                r.amp(3,g)=tmpPar(1)+tmpPar(3);
                % weighted time constant
                r.tDecay(3,g)=(tmpPar(1)*tmpPar(2) + tmpPar(3)*tmpPar(4))/r.amp(3,g);
                % determine width of fitted function:
                % - fit
                curCutoutFit=f(ix3);
                % - indexes to points of fit flanking intersection with zero 
                x1=find(curCutoutFit>r.amp(end,g)*ap.widthFac,1,'last');
                if ~isempty(x1) && x1<=ix(end)
                  x2=x1+1;
                  % corresponding amplitudes
                  y1=curCutoutFit(x1);
                  y2=curCutoutFit(x2);
                  % intersection:
                  r.width(1,g)=x1+(r.amp(end,g)*ap.widthFac-y1)*(x2-x1)/(y2-y1);
                end
                % now embed
                cutoutFit(ix4,g)=curCutoutFit;
                % quality of fit:
                % 1. adjusted r square
                r.qFit(g)=gof.adjrsquare;
                % 2. root mean square error
                r.qFit2(g)=gof.rmse;
                % 3. mean squared error minus mean squared error of noise
                % (as obtained via hipass filtering of traces), normalized
                % to PSC's amplitude
                r.qFit3(g)=(gof.sse - sum(hifiCutout(ix2,g)))/numel(ix)/r.amp(3,g);
              end
              % evaluate fit for time interval stretching up to ten times
              % its decay time into the future and add to dFit
              % - first, index into dFit 
              dFitIx=tsl_pts(g)+r.tPeak(g)+wp.postPeakFitDelay_pts-1;
              tIx=dFitIx:max(dFitIx+10*r.tDecay(end,g),wp.nd);
              dFit(tIx)=dFit(tIx)+f(tIx-tIx(1)+1);
              if gof.adjrsquare>ap.minFitQuality
                % if current fit is good as judged by the adjusted r square
                % criterion embed it in dFitHist. Note that using radj as
                % criterion favors large PSCs.
                dFitHist(tIx)=dFitHist(tIx)+f(tIx-tIx(1)+1);
              end
            catch
              disp(['event # ' int2str(g) ' could not be fitted']);
              r.OKIx(g)=false;
            end
          end
          waitbar(g/wp.coNCol,h);
        end
        delete(h);
        % ** now, mark as not OK those events whose fit error is beyond
        % specified limits
        tmpIx=find(r.qFit<wp.qFitLim | r.qFit3>wp.qFit3Lim);
        r.OKIx(tmpIx)=false;
        disp(['deleting ' int2str(numel(tmpIx)) ' PSCs based on fit error (R squared & normalized error)']);
        % consolidate r
        r=consolidateResults(r);
        % now plot exemplary traces, avoiding extremes and of course
        % cutouts that could not be fitted
        tmpN=16;
        [nix,tmpCoIx]=sort(r.amp(1,:));
        tmpCoIx(isnan(nix))=[];
        tmpCoIx=tmpCoIx(round(linspace(1,numel(tmpCoIx),tmpN+2)));
        tmpCoIx=unique(tmpCoIx(2:end-1));
        plotFitTrace(ap,wp,r,cutout,cutoutFitHist,sp,tmpCoIx);
        if wp.doDumpPSC
          % dump a few juicy pieces on the base workspace:
          assignin('base','d',d);
          assignin('base','dFit',dFit);
          assignin('base','dFitHist',dFitHist);
          assignin('base','cutoutFitHist',cutoutFitHist);
          assignin('base','cutout',cutout);
          assignin('base','cutoutFit',cutoutFit);
          % copy-and-paste line for inspecting pscs and their fits:
          % cix=1; figure(2), cla, plot(cutout(:,cix),'k'), hold on, plot(cutoutFit(:,cix)), plot(cutoutFitHist(:,cix),'c')
        end
        % - finally, update popupmenu UIs
        lst=r2array(ap.nDecayComponent,'parameterList',r);
        set(sp.scatter(1).xAxMenuH,'string',lst,'value',find(strcmp(lst,'time stamp')));
        set(sp.scatter(1).yAxMenuH,'string',lst,'value',find(strcmp(lst,'fit quality:R2adj')));
        % §§ also update overview cutouts plot!
        job(1)={'plotScatter'};
      else
        if isempty(cutout)
          h=warndlg('no data loaded');
        else
          h=warndlg('cutouts must be prepared before fit');
        end
        uiwait(h);
        job(1)=[];
      end        

    case 'plotScatter'
      % create scatter plot of select parameters
      % retrieve data according to current axis selections IN REAL TIME
      % UNITS
      wp.scatterD=r2array(ap.nDecayComponent,...
        [get(sp.scatter(1).xAxMenuH,'value') get(sp.scatter(1).yAxMenuH,'value')],...
        r,tsl,wp.si/1000);
      subplot(sp.scatter(1).axH); 
      % wipe plot and don't forget to delete handles, too
      cla
      wp.curScatCoPh=[];
      ph=plot(wp.scatterD(:,1),wp.scatterD(:,2),'ko');
      % bg color for scatter plot has to be reset
      set(sp.scatter.axH,'color',wp.stdAxCol);      
      % if cutouts had been marked before plot them
      if ~isempty(wp.curScatCoIx)
        switch wp.scatterCursorMode
          case 'probe'
            wp.curScatCoPh=plot(wp.scatterD(wp.curScatCoIx,1),wp.scatterD(wp.curScatCoIx,2),'bo');
          case 'delete'
            wp.curScatCoPh=plot(wp.scatterD(wp.curScatCoIx,1),wp.scatterD(wp.curScatCoIx,2),'ro');
        end
      end
      nicexyax(30);
      set(sp.scatter(1).axH,'buttondownfcn',{@pscfitguifunc,{'reactCursorScatter'}});
      % ** it is important that the data plot have the same callback
      % because otherwise clicks on the points result in naught
      set(ph,'buttondownfcn',{@pscfitguifunc,{'reactCursorScatter'}});
      job(1)=[];
      
    case 'setScatterCursorMode'
      tmp=get(sp.scatter(1).cursorModeMenuH,'string');
      wp.scatterCursorMode=tmp{get(sp.scatter(1).cursorModeMenuH,'value')};
      % reset current cutouts
      wp.curScatCoIx=[];
      wp.curScatCoPh=[];      
      % if user changed scatter mode although no data are to be plotted,
      % don't call plotScatter because it would be pointless
      if ~isempty(wp.scatterD)
        job(1)={'plotScatter'};
      else
        job(1)=[];
      end

    case 'reactCursorScatter'
      % action following mouse click on scatterplot
      if ~isempty(wp.scatterD)
        % start by retrieving axis limits
        xl=get(sp.scatter(1).axH,'xlim');
        yl=get(sp.scatter(1).axH,'ylim');
        % current mouse click
        xy=get(gca,'CurrentPoint');
        % we're only 2D
        xy=xy(1,1:2);
        
        % differentiate according to kind of click:
        switch lower(get(gcf,'SelectionType'))
          case 'extend'
            % shift-click: mark all above current y coordinate
            threshDist= -xy(2);
            dist= -wp.scatterD(:,2);
          case 'alt'
            % ctrl-click: mark all below current y coordinate
            dist=wp.scatterD(:,2);
            threshDist=xy(2);
          otherwise
            % click: swat data points with a circle with a radius of 5% of
            % both axes' extent
            threshDist=.05*diff(xl);
            yScaleFac=diff(xl)/diff(yl);
            % distance of all points from spot clicked on, corrected for
            % differences in axis range
            dist=sqrt((wp.scatterD(:,1)-xy(1)).^2+(yScaleFac*(wp.scatterD(:,2)-xy(2))).^2);
        end
        
        subplot(sp.scatter(1).axH), hold on
        % first, remove previously plotted pts, if any
        if ~isempty(wp.curScatCoPh)
          delete(wp.curScatCoPh);
        end
        % differentiate according to mode
        switch wp.scatterCursorMode
          case 'probe'
            wp.curScatCoIx=find(dist<threshDist);
            % inspecting data in 'probe' mode is tantamount to reverting
            % any choice of cutouts in 'delete' mode, so we have to re-set 
            % r.OKIx each time 
            r.OKIx(1:end)=true;
            tmpPlotString='bo';
          case 'delete'
            % *** philosophy of deletion of events: cumulatively tag
            % (=collect) all clicked-upon events in both wp.curScatCoIx and
            % r.OKIx; do all the important deletions et al in subjob
            % 'deleteEvents'
            wp.curScatCoIx=unique([find(dist<threshDist); wp.curScatCoIx]);
            r.OKIx(wp.curScatCoIx)=false;
            tmpPlotString='ro';            
        end
        % now plot events in color depending on mode
        wp.curScatCoPh=plot(wp.scatterD(wp.curScatCoIx,1),wp.scatterD(wp.curScatCoIx,2),tmpPlotString);
        % ...and also the cutouts (max 20)
        plotFitTrace(ap,wp,r,cutout,cutoutFitHist,sp,wp.curScatCoIx(1:min(20,numel(wp.curScatCoIx))));
        job(1)=[];
      else
        job(1)=[];
      end
      
    case 'deleteEvents'
      if wp.isPrepareFit
        % first of all, retrieve indices of cutouts to be deleted from
        % r.OKIx, NOT from wp.curScatCoIx because wp.curScatCoIx may contain
        % cutouts tagged in 'probe' mode
        delIx=sort(find(~r.OKIx));
        if ~isempty(delIx)
          % - delete ts
          tsl(delIx)=[];
          tsl_pts(delIx)=[];
          %  - take care of entries in r (use r.OKIx)
          r=refreshResults(r);
          % - delete cutouts
          cutout(:,delIx)=[];
          cutoutFilt(:,delIx)=[];
          if ~isempty(cutoutFitHist), cutoutFitHist(:,delIx)=[]; end
          wp.coNCol=size(cutout,2);
          % get rid of markers corresponding to tagged events
          delete(wp.curScatCoPh);
          % - empty wp.curScatCoIx and wp.curScatCoPh because after deletion
          % of the foul events their job is done
          wp.curScatCoIx=[];
          wp.curScatCoPh=[];
          % now plot exemplary traces, avoiding extremes and of course
          % cutouts that could not be fitted
          tmpN=16;
          if wp.isFit
            [nix,tmpCoIx]=sort(r.amp(1,:));
          else
            [nix,tmpCoIx]=sort(r.aPeak(1,:));
          end
          tmpCoIx(isnan(nix))=[];
          tmpCoIx=tmpCoIx(round(linspace(1,numel(tmpCoIx),tmpN+2)));
          tmpCoIx=unique(tmpCoIx(2:end-1));
          if wp.isFit
            plotFitTrace(ap,wp,r,cutout,cutoutFitHist,sp,tmpCoIx);
          else
            plotProcessTrace(wp,r,cutout,cutoutFilt,sp,tmpCoIx);
          end
          disp([int2str(numel(delIx)) ' events were deleted']);
        end
        job(2:end+1)=job;
        job(1:2)={'plotOv','plotScatter'};
      else
        warndlg('saving/deleting events can only be done after pre-processing and/or fitting');
        job(1)=[];
      end

    case 'saveResults'
      % dump ap, ds, and wp into results file so parameters can be
      % retrieved post-analysis
      fitHead.ds=ds;
      fitHead.ap=ap;
      fitHead.wp=wp;
      
      % ** important: remove all graphics handles from variables (wp)
      % because saving graphics handles in Matlab versions 2014b and later
      % implies writing the whole graphics object to file!
      % - first get rid of the major culprit, a struct of handles
      fitHead.wp=rmfield(fitHead.wp,'handles');
      % other than that, fitHead.wp contains only curScatCoPh; other
      % potential handles to (deleted) objects cannot easily be identified
      % (but do no harm, either)
      fitHead.wp=rmfield(fitHead.wp,'curScatCoPh');
      
      % save results under different name, converting to real time units
      fitResult=r;
      fitResult.tPeak=fitResult.tPeak*wp.si/1000;
      fitResult.tRise=fitResult.tRise*wp.si/1000;
      fitResult.tDecay=fitResult.tDecay*wp.si/1000;
      fitResult.width=fitResult.width*wp.si/1000;
      fitResult.xIntvFitEnd=fitResult.xIntvFitEnd*wp.si/1000;
      % also, save tsl because it may differ from original one
      fitResult.tsl=tsl';
      % place results into the 'res' file, not the cutouts file
      save([ds.evtTslFn '.mat'],'fitResult','fitHead','-mat','-append');
      clear fitHead fitResult
      job(1)=[];

    case 'done'
      tmph=findobj('name', 'PSC fit', 'type', 'figure');
      if ~isempty(tmph)
        delete(tmph);
      end
      job(1)=[];
      clear global

    otherwise
      error(['illegal job: ' partJob]);
  end
  done=isempty(job);
end




% ******************** local funcs ***************************************
% ******************** local funcs ***************************************

function adjustXIntPatch(sp,wp)
% adjust patches' x limits: extend their length by si/2 on either side so
% that there are no visual gaps in immediately adjacent intervals
tmpyl=get(sp.cutout.axH,'ylim');
% as the intervals may overlap make the patches transparent
set(sp.cutout.xIntvBaseLinePatchH,'ydata',tmpyl([1 2 2 1])',...
  'xdata',wp.xIntvBaseLine_pts([1 1 2 2])'+[-.5 -.5 .5 .5]', 'facealpha',.8);
set(sp.cutout.xIntvPeakPatchH,'ydata',tmpyl([1 2 2 1])',...
  'xdata',wp.xIntvPeak_pts([1 1 2 2])'+[-.5 -.5 .5 .5]', 'facealpha',.8);
set(sp.cutout.xIntvFitEndMinLineH,'ydata',tmpyl([1 2])',...
  'xdata',wp.xIntvFitEndMin_pts([1 1])'+[.5 .5]');
set(sp.cutout.xIntvFitEndMaxLineH,'ydata',tmpyl([1 2])',...
  'xdata',wp.xIntvFitEndMax_pts([1 1])'+[.5 .5]');

function wp=convertXInt(ds,ap,wp)
% *** compute time variables in pts ***
% we need two different versions of discrete time variables: one for the
% plots, with negative values for pretrigger samples, and another
% for access into variable 'cutouts'.
% So, all *_pts variables are of the former sort, for plotting
% purposes. Field .xOffs_pts represents the offset in points that
% must be SUBTRACTED from these variables in order to access cutouts.
wp.xOffs_pts=cont2discrete(ap.winEvtCutout(1),wp.si/1000,'intv',1);
wp.preTrigDeadT_pts=cont2discrete(ap.preTrigDeadT,wp.si/1000,'intv',1);
wp.xIntvBaseLine_pts=cont2discrete(ap.xIntvBaseLine,wp.si/1000,'intv',1);
wp.xIntvPeak_pts=cont2discrete(ap.xIntvPeak,wp.si/1000,'intv',1);
wp.xIntvFitEndMin_pts=cont2discrete(ap.xIntvFitEndMin,wp.si/1000,'intv',1);
wp.xIntvFitEndMax_pts=cont2discrete(ap.xIntvFitEndMax,wp.si/1000,'intv',1);
wp.postPeakFitDelay_pts=cont2discrete(ap.postPeakFitDelay,wp.si/1000,'intv',1);
% note definition of interval for rise time
wp.xIntvRiseT_pts=[wp.xIntvBaseLine_pts(2)  wp.xIntvPeak_pts(2)];
wp.smoothTSpan_pts=cont2discrete(ap.smoothTSpan,wp.si/1000,'intv',1);
% make sure tspan is odd
if isfinite(wp.smoothTSpan_pts) && ~mod(wp.smoothTSpan_pts,2)
  wp.smoothTSpan_pts=wp.smoothTSpan_pts+1;
end
wp.maxSlowDecay_pts=cont2discrete(ap.maxSlowDecay,wp.si/1000,'intv',0);
% time axis for cutouts in pts for plotting purposes
wp.timeAx_pts=cont2discrete(ap.winEvtCutout,wp.si/1000,'intv',1);
wp.timeAx_pts=wp.timeAx_pts(1):wp.timeAx_pts(2);
if numel(wp.timeAx_pts) ~= wp.coNRow
  error('internal: mismatch between cutout length and x axis');
end

function r=preallocateResults(r,nEvent,nComponent)
% *** fields of are generally 1 by <number of events>. If multiple
% exponentials are fit, row order is
%     fast component | slow component | weighted
% This applies to amplitude and decay time
if nComponent>1
  nRow=nComponent+1;
else
  nRow=1;
end
tmp=nan(1,nEvent);
% index to 'good' events 
r.OKIx=true(1,nEvent);
% base line
r.base=tmp;
% amplitude of peak in filtered event
r.aPeak=tmp;
% time of occurrence of peak in filtered event
r.tPeak=tmp;
% time of occurrence of 10% amplitude in filtered event
r.t10=tmp;
% time of occurrence of 90% amplitude in filtered event
r.t90=tmp;
% rise time computed from the variables above
r.tRise=tmp;
% amplitude of decay component(s)
r.amp=repmat(tmp,nRow,1);
% decay time(s) 
r.tDecay=repmat(tmp,nRow,1);
% width estimated from fit
r.width=tmp;
% right border of fit interval (which may be variable)
r.xIntvFitEnd=tmp;
% quality of fit
r.qFit=tmp;
r.qFit2=tmp;
r.qFit3=tmp;


function r=consolidateResults(r)
% set all entries except OKIx to nan for events which caused trouble during
% either preprocessing or fitting
r.base(~r.OKIx)=nan;
r.aPeak(~r.OKIx)=nan;
r.tPeak(~r.OKIx)=nan;
r.t10(~r.OKIx)=nan;
r.t90(~r.OKIx)=nan;
r.tRise(~r.OKIx)=nan;
r.amp(:,~r.OKIx)=nan;
r.tDecay(:,~r.OKIx)=nan;
r.width(~r.OKIx)=nan;
r.xIntvFitEnd(~r.OKIx)=nan;
r.qFit(~r.OKIx)=nan;
r.qFit2(~r.OKIx)=nan;
r.qFit3(~r.OKIx)=nan;

function r=refreshResults(r)
% delete all entries tagged in .OKIx
r.base(~r.OKIx)=[];
r.aPeak(~r.OKIx)=[];
r.tPeak(~r.OKIx)=[];
r.t10(~r.OKIx)=[];
r.t90(~r.OKIx)=[];
r.tRise(~r.OKIx)=[];
r.amp(:,~r.OKIx)=[];
r.tDecay(:,~r.OKIx)=[];
r.width(~r.OKIx)=[];
r.xIntvFitEnd(~r.OKIx)=[];
r.qFit(~r.OKIx)=[];
r.qFit2(~r.OKIx)=[];
r.qFit3(~r.OKIx)=[];
% last
r.OKIx(~r.OKIx)=[];

function [wp,r]=discardResults(wp,sp)
wp.isPrepareFit=false;
wp.isFit=false;
r=[];
% wipe plots
cla(sp.selectCutout.axH);
cla([sp.scatter.axH])
wp.scatterD=[];
wp.curScatCoPh=[];
wp.curScatCoIx=[];
set(sp.scatter(1).xAxMenuH,'string','x axis selection','value',1);
set(sp.scatter(1).yAxMenuH,'string','y axis selection','value',1);
drawnow;

function plotProcessTrace(wp,r,cutout,cutoutFilt,sp,coIx)
% plot a selection of cutouts, filtered versions and base line, peak and
% 10 and 90% amplitude points, if available
if ~isempty(coIx)
  % - portion of cutout (note that 1st point is identical to ix as
  % defined for detection of rise time above, facilitating plotting)
  ix=(wp.xIntvRiseT_pts(1):ceil(mean([wp.xIntvPeak_pts(2) wp.xIntvFitEndMin_pts])))-wp.xOffs_pts;
  tmpN=numel(coIx);
  subplot(sp.selectCutout.axH), cla
  hold on
  % unfiltered ones first
  [~,dy]=pllplot(cutout(ix,coIx),'spacing','percentile','dy',10,'noplot',1);
  [ylim,dy,~,ph]=pllplot(cutout(ix,coIx),'spacing','fixed',...
    'dy',dy);
  set(ph,'linewidth',2.5);
  [~,dy,~,ph]=pllplot(cutoutFilt(ix,coIx),'spacing','fixed',...
    'dy',dy,'ylim',ylim);
  set(ph,'color',wp.FitEventCol,'linewidth',1.5);
  % plot detected peaks
  tmpAmp=r.aPeak(coIx)+cumsum([0 dy]);
  ph=plot(r.tPeak(coIx)-wp.xIntvRiseT_pts(1)+1,tmpAmp,'k^');
  set(ph,'markerfacecolor',wp.xIntvPeakCol,'color',wp.xIntvPeakCol,'markersize',8);
  % ...base line
  tmpAmp=cumsum([0 dy]);
  ph=plot(wp.xIntvBaseLine_pts(2)-wp.xIntvRiseT_pts(1)+1,tmpAmp,'kv');
  set(ph,'markerfacecolor',wp.xIntvBaseLineCol,'markersize',6);
  % as well as 10 and 90% amplitude pts
  m1=nan(1,tmpN);
  m2=m1;
  for g=1:tmpN
    tmpX=floor(r.t10(coIx(g)))+[0 1]-wp.xOffs_pts;
    tmpY=cutoutFilt(tmpX,coIx(g));
    m1(g)=interp1(tmpX,tmpY,r.t10(coIx(g))-wp.xOffs_pts);
    tmpX=floor(r.t90(coIx(g)))+[0 1]-wp.xOffs_pts;
    tmpY=cutoutFilt(tmpX,coIx(g));
    m2(g)=interp1(tmpX,tmpY,r.t90(coIx(g))-wp.xOffs_pts);
  end
  plot(r.t10(coIx)-wp.xIntvRiseT_pts(1)+1,m1+cumsum([0 dy]),'w+');
  plot(r.t90(coIx)-wp.xIntvRiseT_pts(1)+1,m2+cumsum([0 dy]),'w+');
  % if there are three or fewer cutouts to plot squeeze y axis
  if numel(coIx)<=3
    nicexyax(2);
  else
    nicexyax;
  end
  hold off
  axis on
  set(gca,'xtick',[]);
  % set proper background color
  set(sp.selectCutout.axH,'color',wp.stdAxCol);
end

function plotFitTrace(ap,wp,r,cutout,cutoutFitHist,sp,coIx)
% plot a selection of cutouts in raw and fitted versions
% ** has to be accomplished in a loop because start and stop times of fits
% vary from event to event
subplot(sp.selectCutout.axH), cla
set(sp.selectCutout.axH,'color',wp.stdAxCol);
hold on
yOffs=0;
for g=coIx(:)'
  % offset for current cutout
  yOffs=yOffs-r.aPeak(g)*1.1;
  % time axis for fit: unit is still ticks, but we start at zero
  ix=(0:r.xIntvFitEnd(g)-r.tPeak(g)-wp.postPeakFitDelay_pts)';
  % the corresponding index into cutout
  ix2=(r.tPeak(g)+wp.postPeakFitDelay_pts:r.xIntvFitEnd(g))'-wp.xOffs_pts;
  % whole cutout first
  ph=plot(wp.timeAx_pts,cutout(:,g)+yOffs,'k');
  % fit
  % ph=plot(wp.timeAx_pts,cutoutFitHist(:,g)+yOffs,'k');
  switch ap.nDecayComponent
    case 1
      ph=plot(wp.timeAx_pts(ix2),r.amp(1,g)*exp(-ix/r.tDecay(1,g))+yOffs+...
        cutoutFitHist(ix2,g)+r.base(g));
    case 2
      ph=plot(wp.timeAx_pts(ix2),r.amp(1,g)*exp(-ix/r.tDecay(1,g))+...
        r.amp(2,g)*exp(-ix/r.tDecay(2,g))+yOffs+cutoutFitHist(ix2,g)+r.base(g));
  end
  set(ph,'linewidth',1.5,'color',wp.FitEventCol);
end
% if there are three or fewer cutouts to plot squeeze y axis 
if numel(coIx)<=3
  niceyax(2);
else
  niceyax;
end
set(sp.selectCutout.axH,'xlim',wp.timeAx_pts([1 end]));
hold off



function out=r2array(nComponent,item,r,varargin)
% this function serves a dual purpose:
% i. it may output a list of all results parameters (that can be used e.g.
% in listbox UIs)
% ii. it may output values of those parameters (that is, it maps strings to
% fields of results variable r)
% Details:
% - input 'item' must be a char array 'parameterList' or a cell array of
% chars
% - input var 'item' may also be an integer, namely, the index into cell
% array 'list' below. 
% - varargin{1}=tsl
% - varargin{2}: sampling interval in arbitrary units. If specified, 
% parameters tPeak, tRise and tDecay will be given in these time units
% ** when introducing new parameters with time as unit make sure conversion
% is implemented here AND in the 'save results' job 
if nargin==5
  tConvFac=varargin{2};
else
  tConvFac=1;
end

list1={'tau','amplitude'};
list2={...
'tau (fast)',...
'amplitude (fast)',...
'tau (slow)',...
'amplitude (slow)',...
'tau (weigh)',...
'amplitude (weigh)'...
};
list_common={...
'width',...
'10-90% rise time',...
'fit quality:R2adj',...
'fit error:rmse',...
'fit error:MSE(non-noise)/amp',...
'fit interval',...
'baseline',...
'time stamp'...
};

switch nComponent
  case 1
    list=cat(2,list1,list_common);
  case 2
    list=cat(2,list2,list_common);
end

if ischar(item) && strcmp(item,'parameterList')
  out=list;
else
  out=nan(numel(r.base),numel(item));
  % if 'item' is a numeric array convert into cell array
  if isnumeric(item)
    item=list(item);
  end
  % now 'translate' strings into results fields
  for g=1:numel(item)
    switch item{g}
      % one component:
      case 'tau'
        out(:,g)=r.tDecay'*tConvFac;
      case 'amplitude'
        out(:,g)=r.amp';
        % two components:
      case 'tau (fast)'
        out(:,g)=r.tDecay(1,:)'*tConvFac;
      case 'amplitude (fast)'
        out(:,g)=r.amp(1,:)';
      case 'tau (slow)'
        out(:,g)=r.tDecay(2,:)'*tConvFac;
      case 'amplitude (slow)'
        out(:,g)=r.amp(2,:)';
      case 'tau (weigh)'
        out(:,g)=r.tDecay(3,:)'*tConvFac;
      case 'amplitude (weigh)'
        out(:,g)=r.amp(3,:)';
        % common parameters:
      case 'width'
        out(:,g)=r.width(1,:)'*tConvFac;
      case '10-90% rise time'
        out(:,g)=r.tRise'*tConvFac;
      case 'fit quality:R2adj'
        out(:,g)=r.qFit';
      case 'fit error:rmse'
        out(:,g)=r.qFit2';
      case 'fit error:MSE(non-noise)/amp'
        out(:,g)=r.qFit3';
      case 'fit interval'
        out(:,g)=r.xIntvFitEnd'*tConvFac; 
      case 'baseline'
        out(:,g)=r.base';
      case 'time stamp'
        out(:,g)=varargin{1};
      otherwise
        error(['internal:illegal item for ' mfilename]);
    end
  end
end
