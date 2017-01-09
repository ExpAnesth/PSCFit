function pscfitguifunc(src,eventdata,job,varargin)
% ** function pscfitguifunc(src,eventdata,job,varargin)
% Collection of callback routines for pscfitgui.m
%                         >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% job              cell array of char    jobs to accomplish
% sp               struct                handles to subplots of main
%                                         figure


% We need persistent variables:
%   tsl = cell array of time stamp lists of single 'events'
%   cutout = cutouts of events
%   fCutout = filtered version (solely for peak detection and plotting
%   purposes)
%   r = results structure
% Strucures:
%   ds='data set' listing properties of current file
%   ap='analysis parameters' describing details of current analyis
%   wp='working parameters' (like colors)
%   sp=subplot handles
persistent ds ap wp sp tsl tsl_pts cutout fCutout head r


% to do:
% - conversion to cont time
% - possibly make results var r an array with fixed column indices (baseIx,...) ?


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
      ds.resFn='lastFile_res';
      % name of file holding evtCutout
      ds.evtCutoutFn='lastFile_evtCutout';
      ds.dataFnCore='lastFile';
      ds.dataPath='d:\_data\_IPSCFit\';
      % cutout window (ms)
      ds.winEvtCutout=[nan nan];
      % +1 for positive-going events, -1 for negative-going events
      ds.polarity=nan;
      % ds.fileInfo=[];

      % -------------------------------------
      % ----- set up ap (analysis parameters)
      % -------------------------------------
      % ~~~~~~~ file names
      % ~~~~~~~ standard parameters section
      % - number of exponential components to fit to current decay
      ap.nDecayComponent=1;
      % - pre-trigger dead time (ms) (the posttrigger dead time is
      % equivalent to xIntvFitEndMin below
      ap.preTrigDeadT=nan;
      % - interval in which to compute base line (ms)
      ap.xIntvBaseLine=[nan nan];
      % - interval in which to seek peak (ms)
      ap.xIntvPeak=[nan nan];
      % - minimally acceptable value of right border of interval in which to fit (ms)
      ap.xIntvFitEndMin=nan;
      % - maximal value of right border of interval in which to fit (ms)
      ap.xIntvFitEndMax=nan;
      % - fitting interval starts postPeakFitDelay ms after peak
      ap.postPeakFitDelay=nan;
      % - was average of cutouts processed?
      ap.isAverageCutout=0;
      % - comment on set of parameters
      ap.parComment='a little space for notes';
      % ~~~~~~~ advanced parameters section
      % - smoothing of data: span of filter (ms) (span=number of neighbours
      % of a data point that shall make it into the smoothing process)
      ap.smoothTSpan=nan;
      % - applies to data smoothing with Savitzky-Golay filters: polynomial
      % order
      ap.smoothPolyOrder=5;

      % --------------------------------------
      % ----- set up wp ('working' parameters)
      % --------------------------------------
      % ~~~~~~~ display options & matlab version section
      % which version of matlab?
      wp.mver=ver;
      wp.mver=str2double(wp.mver(strmatch('matlab',lower({wp.mver.Name}),'exact')).Version);
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
      wp.resFnString='';
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
      % --------------------------------------
      % ----- set up results struct
      % --------------------------------------
      r=preallocateResults(r,1,1);

      % --------------------------------------
      % ----- initialize subplots
      % --------------------------------------
      % -- cutouts (all)
      % bg and cutout colors 
      set(sp.cutout.axH,'color',wp.stdAxCol,'colororder',coma('bluered','ncols',50));
      subplot(sp.cutout.axH), hold on,
      % plot the patches representing x intervals
      tmpyl=get(sp.cutout.axH,'ylim');
      sp.cutout.xIntvBaseLinePatchH=patch(wp.xIntvBaseLine_pts([1 1 2 2])',tmpyl([1 2 2 1])',wp.xIntvBaseLineCol);
      sp.cutout.xIntvPeakPatchH=patch(wp.xIntvPeak_pts([1 1 2 2])',tmpyl([1 2 2 1])',wp.xIntvPeakCol);
      % the lines representing minimal and maximal right borders of fit interval
      sp.cutout.xIntvFitEndMinLineH=line([nan nan]',tmpyl([1 2])','color',wp.xIntvFitEndMinCol,'LineWidth',1.5);
      sp.cutout.xIntvFitEndMaxLineH=line([nan nan]',tmpyl([1 2])','color',wp.xIntvFitEndMaxCol,'LineWidth',1.5);
      % handle to cutouts
      sp.cutout.ph=plot([nan]);
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
      uicFn=fieldnames(wp.handles);
      % apFn contains names of all fields that can be set by uicontrols
      apFn=intersect(fieldnames(ap),uicFn);
      % same for wp
      wpFn=intersect(fieldnames(wp),uicFn);
      % by default look for files in \PSCFit\parms
      pfDir=mfilename('fullpath');
      pfDir=pfDir(1:max(strfind(pfDir,'\'))-1);
      w=what;
      cd([pfDir '\parms']);
      [tmpOptFn,tmpOptPath] = uigetfile('*.mat','pick parameter file');
      if ischar(tmpOptFn) && ischar(tmpOptPath)
        load([tmpOptPath tmpOptFn]);
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
      % cd back to original dir
      cd(w.path);
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
        structIx=[~isempty(strmatch(uicFn{g},apFn,'exact')),...
          ~isempty(strmatch(uicFn{g},wpFn,'exact'))];
        if length(find(structIx))==1
          eval(['cType=get(wp.handles.' uicFn{g} ',''style'');']);
          switch cType
            case 'edit'
              % depending on the type of the field...
              switch uicFn{g}
                case {'parComment','resFnString'}
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
        structIx=[~isempty(strmatch(uicFn{g},apFn,'exact')),...
          ~isempty(strmatch(uicFn{g},wpFn,'exact'))];
        if length(find(structIx))==1
          eval(['cType=get(wp.handles.' uicFn{g} ',''style'');']);
          switch cType
            case 'edit'
              % depending on the type of the field...
              switch uicFn{g}
                case {'parComment','resFnString'}
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
      % fields of ap and wp. For example, resFn is a field of ap that
      % cannot be set via a uicontrol but depends on the specific file that
      % is currently analyzed. Such file-specific information must not be
      % saved. Instead, only those fields of ap and wp represented by
      % uicontrols shall be saved. This is done in the lines below
      uicFn=fieldnames(wp.handles);
      % ap_uic contains all fields that can be set by uicontrols
      ap_uic=rmfield(ap,setdiff(fieldnames(ap),uicFn));
      % same for wp
      wp_uic=rmfield(wp,setdiff(fieldnames(wp),uicFn));
      % by default dump files in \PSCFit\parms
      pfDir=mfilename('fullpath');
      pfDir=pfDir(1:max(strfind(pfDir,'\'))-1);
      w=what;
      cd([pfDir '\parms']);
      [tmpDataFn,tmpDataPath] = uiputfile('*.mat');
      if ischar(tmpDataFn) && ischar(tmpDataPath)
        save([tmpDataPath tmpDataFn],'ap_uic','wp_uic');
      end
      % cd back to original dir
      cd(w.path);
      % delete vars
      clear w pfDir ap_uic wp_uic
      job(1)=[];

    case 'digestParameters'
      % this part job is run whenever parameters were read from gui
      job(1)=[];
      isAllParameterOK=true;
      % if any parameter changes after events have been fit previous
      % results must be discarded 
      doDigest=true;
      butt='Yes';
      % §§§ problem with approach below: the changed parameter value(s)
      % cannot be reverted in the current implementation in which old
      % values are not stored. For now, kill results without making a big fuss.
%       if wp.isFit || wp.isPrepareFit
%         butt=questdlg('Changing parameters will discard results from previous fit - continue?','question','Yes');
%         if strcmp(butt,'Yes')
%           [wp,r]=discardResults(wp,r,sp);
%         else
%           doDigest=false;
%         end
%       end

      [wp,r]=discardResults(wp,r,sp);
      if doDigest
        disp('** processing & checking options..');
        % ----- checks of parameters:
        % --- number of decay components
        if isempty(intersect(ap.nDecayComponent, [1 2]))
          warndlg('number of decay components must be 1 or 2 - setting to 1');
          isAllParameterOK=false;;
          ap.nDecayComponent=1;
        end
        % --- x intervals
        % - dead time
        if length(ap.preTrigDeadT)~=1 || ~isfinite(ap.preTrigDeadT)
          warndlg('pre-trigger dead time must be one finite, non-nan value')
          isAllParameterOK=false;
          ap.preTrigDeadT=nan;
        end
        % - base line
        if length(ap.xIntvBaseLine)~=2 || any(~isfinite(ap.xIntvBaseLine))
          warndlg('base line interval must contain two finite, non-nan values')
          isAllParameterOK=false;;
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
        % --- smoothing parameters
        if ap.smoothTSpan<0 || ~isfinite(ap.smoothTSpan)
          warndlg('smoothing time span must be finite and >0 - setting to 0.1');
          isAllParameterOK=false;
          ap.smoothTSpan=.1;
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

      if ~isAllParameterOK || ~strcmp(butt,'Yes')
        job(1)={'writeParameters2Gui'};
      end

    case 'readData'
      % ******************************************************************
      % delete/reset results/plots/variables from last file/channel
      cutout=[];
      fCutout=[];
      wp.si=nan;
      tsl=[];
      ds.winEvtCutout=[nan nan];
      wp.timeAx_pts=nan;
      ds.polarity=nan;
      ap.isAverageCutout=0;
      % this takes care of all fit results and plots except sp.cutout.axH 
      [wp,r]=discardResults(wp,r,sp);
      if any(ishandle(sp.cutout.ph))
        delete(sp.cutout.ph);
        sp.cutout.ph=nan;
      end
      set(sp.cutout.th,'string','no data loaded');

      [tmpDataFn,tmpDataPath] = uigetfile('*.mat','pick *evtCutout* file',ds.dataPath);
      % ** note that only *evtCutout* files will be accepted
      tmpCond=[ischar(tmpDataFn) ischar(tmpDataPath)  ~isempty(strfind(tmpDataFn,'evtCutout'))];
      if isequal(tmpCond,[true true false])
        warndlg('please pick a *evtCutout*.mat file');
      end
      if all(tmpCond)
        % names of .mat files:
        tmpix=strfind(tmpDataFn,'_evtCutout');
        % - core file name
        ds.dataFnCore=tmpDataFn(1:tmpix(end)-1);
        % - path
        ds.dataPath=tmpDataPath;
        % - full name of file containing the cutouts (the one picked above)
        ds.evtCutoutFn=[tmpDataPath tmpDataFn];
        ds.resFn=[tmpDataPath ds.dataFnCore '_res'];
        % load cutouts and head
        load(ds.evtCutoutFn);
        % load time stamps
        load(ds.resFn,'evt');
        if numel(evtCutout)>1
          error('evtCutout contains more than one set of cutouts');
        end
        % *** copy into local variables
        cutout=evtCutout{1};
        tsl=evt.tsl{1};
        clear evtCutout evt
        % size of things
        [wp.coNRow,wp.coNCol]=size(cutout);
        % the original cutout window (ms)
        ds.winEvtCutout=head.ap.winEvtCutout;
        % also, copy sampling interval to wp
        wp.si=head.ds.fileInfo.si;
        % now convert time intervals to pts
        wp=convertXInt(ds,ap,wp);
        % same for time stamps
        tsl_pts=cont2discrete(tsl,wp.si/1000,'intv',0);
        % next job: plot data
        job(1)={'plotOv'};
        clear tmp*
      else
        job(1)=[];
      end

    case 'plotOv'
      if ~isempty(cutout)
        % plot specified number of cutouts
        tmpv1=min(wp.maxNPlotCutout,wp.coNCol);
        tmpix=unique(ceil((1:tmpv1)/tmpv1*wp.coNCol));
        subplot(sp.cutout.axH)
        % - set x and y values of cutout data to new values
        set(sp.cutout.axH,...
          'xlim',wp.timeAx_pts([1 end]),...
          'ylim',[min(min(cutout(:,tmpix)))  max(max(cutout(:,tmpix)))],'colororder',coma('bluered','ncols',50));
        sp.cutout.ph=plot(wp.timeAx_pts,cutout(:,tmpix));
        % set real time units as labels (such that zero and integer values
        % are shown. This can be achieved easily by shifting x ticks as
        % defined by wp.timeAx_pts one tick to the right)
        tmpxtick=cont2discrete(5,wp.si/1000,'intv',1);
        tmpxtick=[fliplr([1:-tmpxtick:wp.timeAx_pts(1)]) 1+tmpxtick:tmpxtick:wp.timeAx_pts(end)];
        set(sp.cutout.axH,'xtick',tmpxtick);
        set(sp.cutout.axH,'xticklabel',discrete2cont(tmpxtick,wp.si/1000,'intv',0));
        adjustXIntPatch(sp,wp);
        % information on file
        if ishandle(sp.cutout.th)
          delete(sp.cutout.th);
        end
        sp.cutout.th=smarttext([ds.dataFnCore ', ' int2str(tmpv1) ' of ' int2str(wp.coNCol) ' events'],...
          .99,.03,'color','k','fontsize',10,'fontweight','bold','interpreter','none');
        if ap.isAverageCutout
          set(sp.cutout.th,'string',[ds.dataFnCore ', AVERAGE']);
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
          fCutout=[];
          % this takes care of all fit results and plots except sp.cutout.axH
          [wp,r]=discardResults(wp,r,sp);
          if any(ishandle(sp.cutout.ph))
            delete(sp.cutout.ph);
            sp.cutout.ph=nan;
          end
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

    case  'prepareFit'
      if ~isempty(cutout)
        wp.isPrepareFit=true;
        r=preallocateResults(r,wp.coNCol,ap.nDecayComponent);
        % wipe plots
        cla(sp.selectCutout.axH);
        cla([sp.scatter.axH])
        wp.scatterD=[];
        wp.curScatCoPh=[];
        wp.curScatCoIx=[];
        % i) base line 
        ix=(wp.xIntvBaseLine_pts(1):wp.xIntvBaseLine_pts(2))-wp.xOffs_pts;
        r.base=mean(cutout(ix,:));
        % keep values in .base but subtract from traces
        cutout=cutout-repmat(r.base,wp.coNRow,1);
        % ii) polarity
        % if average of points in peak search interval positive we're
        % dealing with positive-going events and vice versa
        ix=(wp.xIntvPeak_pts(1):wp.xIntvPeak_pts(2))-wp.xOffs_pts;
        ds.polarity=sign(median(median(cutout(ix,:))));
        % **************************************************************
        % invert negative-going events to streamline all following 
        % operations 
        % **************************************************************
        if ds.polarity<0
          cutout=cutout*-1;
        end

        % iii) peaks, following Savitzky-Golay filtering
        % - keep current ix 
        fCutout=sgolayfilt(cutout,ap.smoothPolyOrder,wp.smoothTSpan_pts);
        % ** note that tPeak is in pts, corrected for xIntvPeak interval
        % such that continuous time zero corresponds to discrete time 1, as
        % all _pts parameters
        evr=evdeal(fCutout(ix,:),'idx',{'minmaxpeak'});
        r.aPeak=evr.maxPeak;
        r.tPeak=evr.maxPeakT+ix(1)-1+wp.xOffs_pts;
        
        % iv) 10-90% rise time (determine from filtered version, too)
        ix=(wp.xIntvRiseT_pts(1):wp.xIntvRiseT_pts(2))-wp.xOffs_pts;
        [m,i1]=max(fCutout(ix,:)>=repmat(.1*r.aPeak,numel(ix),1));
        [m,i2]=max(fCutout(ix,:)>=repmat(.9*r.aPeak,numel(ix),1));
        r.tRise=i2-i1;
        % adjust time frame as above
        i1=i1+ix(1)-1+wp.xOffs_pts;
        i2=i2+ix(1)-1+wp.xOffs_pts;
        
        % v) tag all events which are problematic:
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
        
        % vi) plot a selection of tmpN cutouts in raw and filtered versions
        % so we can judge whether filters pars are appropriate. Take events
        % from whole pool of cutouts, not just the ones deemed OK.
        % - portion of cutout (note that 1st point is identical to ix as
        % defined for detection of rise time above, facilitating plotting)
        ix=(wp.xIntvRiseT_pts(1):ceil(mean([wp.xIntvPeak_pts(2) wp.xIntvFitEndMin_pts])))-wp.xOffs_pts;
        % - pick cutouts representing range of amplitudes, avoiding
        % extremes
        tmpN=20;
        [nix,tmpCoIx]=sort(abs(r.aPeak));
        tmpCoIx=tmpCoIx(round(linspace(1,wp.coNCol,tmpN+2)));
        tmpCoIx=unique(tmpCoIx(2:end-1));
        % there may be fewer than 20 coututs
        tmpN=numel(tmpCoIx);
        subplot(sp.selectCutout.axH), cla
        hold on
        % unfiltered ones first
        [ylim,dy]=pllplot(cutout(ix,tmpCoIx),'noplot',1);
        [ylim,dy,yscaleFac,ph]=pllplot(cutout(ix,tmpCoIx),'spacing','fixed',...
          'dy',dy);
        set(ph,'linewidth',2.5);
        [ylim,dy,yscaleFac,ph]=pllplot(fCutout(ix,tmpCoIx),'spacing','fixed',...
          'dy',dy,'ylim',ylim);
        set(ph,'color',wp.FitEventCol,'linewidth',1.5);
        % plot detected peaks
        tmpAmp=evr.maxPeak(tmpCoIx)+cumsum([0 dy]);
        ph=plot(r.tPeak(tmpCoIx)-wp.xIntvRiseT_pts(1)+1,tmpAmp,'k^');
        set(ph,'markerfacecolor',wp.xIntvPeakCol,'color',wp.xIntvPeakCol,'markersize',8);
        % ...base line
        tmpAmp=cumsum([0 dy]);
        ph=plot(wp.xIntvBaseLine_pts(2)-wp.xIntvRiseT_pts(1)+1,tmpAmp,'kv');
        set(ph,'markerfacecolor',wp.xIntvBaseLineCol,'markersize',6);
        % as well as 10 and 90% amplitude points
        for g=1:tmpN
          m1(g)=fCutout(i1(tmpCoIx(g))-wp.xOffs_pts,tmpCoIx(g));
          m2(g)=fCutout(i2(tmpCoIx(g))-wp.xOffs_pts,tmpCoIx(g));
        end
        ph=plot(i1(tmpCoIx)-wp.xIntvRiseT_pts(1)+1,m1+cumsum([0 dy]),'w+');
        % set(ph,'markerfacecolor',mean([wp.xIntvBaseLineCol;wp.xIntvPeakCol]));
        ph=plot(i2(tmpCoIx)-wp.xIntvRiseT_pts(1)+1,m2+cumsum([0 dy]),'w+');
        % set(ph,'markerfacecolor',mean([wp.xIntvBaseLineCol;wp.xIntvPeakCol]));
        nicexax;
        hold off
        axis on
        set(gca,'xtick',[]);
        % set proper background color
        set(sp.selectCutout.axH,'color',wp.stdAxCol);
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
            
          case 2
            ft=fittype('a1*exp(-t/tau1) + a2*exp(-t/tau2)',...
              'independent',{'t'},...
              'coefficients',{'a1','tau1','a2','tau2'});
            % upper limits: make plausibility assumptions
            % - either amplitude cannot be more than x*maximal peak found
            % - decay time is extremely unlikely to be more than
            % 2*posttrigger cutout length
            u=[max(r.aPeak)*2  wp.timeAx_pts(end)*2  max(r.aPeak)*2  wp.timeAx_pts(end)*2];
            % lower bounds: 
            l=[0 0 0 0];
        end
        fo=fitoptions('method','NonlinearLeastSquares',...
          'upper',u,...
          'lower',l,...
          'TolX',.001,...
          'display', 'off'...
          );
        % compute right border of fit interval as described below and store
        % this information in field of r
        r.xIntvFitEnd=min(wp.xIntvFitEndMax_pts,[diff(tsl_pts)' inf]);
        h=waitbar(0,'Fitting...');
        for g=1:wp.coNCol
          if r.OKIx(g)
            % first, the interval to be fitted must be determined:
            % - the first point to be fitted is wp.postPeakFitDelay_pts
            % after the peak
            % - the last point is, ideally, xIntvFitEndMax_pts
            % - if however, there is a following event in this interval,
            % move the interval border back towards zero, but only up to
            % xIntvFitEndMin_pts (all events with a closer event have been
            % kicked out in job preparefit)
            % time axis for fit: unit is still ticks, but we start at zero
            ix=(0:r.xIntvFitEnd(g)-r.tPeak(g)-wp.postPeakFitDelay_pts)';
            % the corresponding index into cutout
            ix2=(r.tPeak(g)+wp.postPeakFitDelay_pts:r.xIntvFitEnd(g))'-wp.xOffs_pts;
            try
              if ap.nDecayComponent==1
                % start values for fit: deduce from parameters computed thus far
                set(fo,'Startpoint',[r.aPeak(g) wp.xIntvFitEndMin_pts]);
                [f,gof]=fit(ix,cutout(ix2,g),ft,fo);
                tmpPar=coeffvalues(f);
                % §§§ confints, predint?
                % amplitude of decay component
                r.amp(1,g)=tmpPar(1);
                % decay time(s)
                r.tDecay(1,g)=tmpPar(2);
                % quality of fit
                r.qFit(g)=gof.adjrsquare;
              else
                % by default, assume that fast component has higher amplitude
                set(fo,'Startpoint',[r.aPeak(g)*.7  wp.xIntvFitEndMin_pts/2   r.aPeak(g)*.3  wp.xIntvFitEndMin_pts*2]);
                [f,gof]=fit(ix,cutout(ix2,g),ft,fo);
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
                % quality of fit
                r.qFit(g)=gof.adjrsquare;
                % summed amplitudes
                r.amp(3,g)=tmpPar(1)+tmpPar(3);
                % weighted time constant
                r.tDecay(3,g)=(tmpPar(1)*tmpPar(2) + tmpPar(3)*tmpPar(4))/r.amp(3,g);
              end
            catch
              disp(['event # ' int2str(g) ' could not be fitted']);
              r.OKIx(g)=false;
            end
          end
          waitbar(g/wp.coNCol,h);
        end
        delete(h)
        % consolidate r
        r=consolidateResults(r);
        % now plot exemplary traces, avoiding extremes and of course
        % cutouts that could not be fitted
        tmpN=20;
        [nix,tmpCoIx]=sort(r.amp(1,:));
        tmpCoIx(isnan(nix))=[];
        tmpCoIx=tmpCoIx(round(linspace(1,numel(tmpCoIx),tmpN+2)));
        tmpCoIx=unique(tmpCoIx(2:end-1));
        plotFitTrace(ap,wp,r,cutout,sp,tmpCoIx);
        % - finally, update popupmenu UIs
        lst=r2array(ap.nDecayComponent,'parameterList',r);
        set(sp.scatter(1).xAxMenuH,'string',lst,'value',1);
        set(sp.scatter(1).yAxMenuH,'string',lst,'value',2);
      else
        if isempty(cutout)
          h=warndlg('no data loaded');
        else
          h=warndlg('cutouts must be prepared before fit');
        end
        uiwait(h);
      end        
      job(1)={'plotScatter'};

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
            % §§§
        end
      end
      axis tight
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
      job(1)=[];

    case 'reactCursorScatter'
      % action following mouse click on scatterplot
      if ~isempty(wp.scatterD)
        % start by retrieving axis limits
        xl=get(sp.scatter(1).axH,'xlim');
        yl=get(sp.scatter(1).axH,'ylim');
        % swat data points with a circle with a radius of 5% of both axes' extent
        threshDist=.05*diff(xl);
        yScaleFac=diff(xl)/diff(yl);
        % current mouse click
        xy=get(gca,'CurrentPoint');
        % we're only 2D
        xy=xy(1,1:2);
        % distance of all points from spot clicked on, corrected for
        % differences in axis range
        dist=sqrt((wp.scatterD(:,1)-xy(1)).^2+(yScaleFac*(wp.scatterD(:,2)-xy(2))).^2);
        tmpCoIx=find(dist<threshDist);

        subplot(sp.scatter(1).axH), hold on
        % accumulate in kill mode but not in probe mode
        switch wp.scatterCursorMode
          case 'probe'
            wp.curScatCoIx=tmpCoIx;
            % first, remove previously plotted pts, if any
            if ~isempty(wp.curScatCoPh) 
              delete(wp.curScatCoPh);
            end
            % now plot them in different color
            wp.curScatCoPh=plot(wp.scatterD(wp.curScatCoIx,1),wp.scatterD(wp.curScatCoIx,2),'bo');
            % ...and also the cutouts (max 20)
            plotFitTrace(ap,wp,r,cutout,sp,wp.curScatCoIx(1:min(20,numel(wp.curScatCoIx))));
          case 'delete'
            % §§ to come
        end
      end
      job(1)=[];      
      
    case 'saveResults'
      % dump ap, ds, and wp into results file so parameters can be
      % retrieved post-analysis
      fitHead.ds=ds;
      fitHead.ap=ap;
      fitHead.wp=wp;
      % save results under different name, converting to real time units
      fitResult=r;
      fitResult.tPeak=fitResult.tPeak*wp.si/1000;
      fitResult.tRise=fitResult.tRise*wp.si/1000;
      fitResult.tDecay=fitResult.tDecay*wp.si/1000;
      fitResult.xIntvFitEnd=fitResult.xIntvFitEnd*wp.si/1000;
      % place results into the 'res' file, not the cutouts file
      save([ds.resFn '.mat'],'fitResult','fitHead','-mat','-append');
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
wp.xOffs_pts=cont2discrete(ds.winEvtCutout(1),wp.si/1000,'intv',1);
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
% time axis for cutouts in pts for plotting purposes
wp.timeAx_pts=cont2discrete(ds.winEvtCutout,wp.si/1000,'intv',1);
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
tmp=repmat(nan,1,nEvent);
% index to 'good' events 
r.OKIx=repmat(true,1,nEvent);
% base line
r.base=tmp;
% amplitude of peak in filtered event
r.aPeak=tmp;
% time of occurrence of peak in filtered event
r.tPeak=tmp;
% rise time
r.tRise=tmp;
% amplitude of decay component(s)
r.amp=repmat(tmp,nRow,1);
% decay time(s) 
r.tDecay=repmat(tmp,nRow,1);
% right border of fit interval (which may be variable)
r.xIntvFitEnd=tmp;
% quality of fit
r.qFit=tmp;

function r=consolidateResults(r)
% set all entries to nan for events which caused trouble during either
% preprocessing or fitting
r.base(~r.OKIx)=nan;
r.aPeak(~r.OKIx)=nan;
r.tPeak(~r.OKIx)=nan;
r.tRise(~r.OKIx)=nan;
r.amp(:,~r.OKIx)=nan;
r.tDecay(:,~r.OKIx)=nan;
r.xIntvFitEnd(~r.OKIx)=nan;
r.qFit(~r.OKIx)=nan;


function [wp,r]=discardResults(wp,r,sp)
wp.isPrepareFit=false;
wp.isFit=false;
r=[];
% wipe plots
cla(sp.selectCutout.axH);
cla([sp.scatter.axH])
wp.scatterD=[];
wp.curScatCoPh=[];
wp.curScatCoIx=[];
drawnow;


function plotFitTrace(ap,wp,r,cutout,sp,coIx)
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
  switch ap.nDecayComponent
    case 1
      ph=plot(wp.timeAx_pts(ix2),r.amp(1,g)*exp(-ix/r.tDecay(1,g))+yOffs);
    case 2
      ph=plot(wp.timeAx_pts(ix2),r.amp(1,g)*exp(-ix/r.tDecay(1,g)) + ...
        r.amp(2,g)*exp(-ix/r.tDecay(2,g))+yOffs);
  end
  set(ph,'linewidth',1.5,'color',wp.FitEventCol);
end
% if there are three or fewer cutouts to plot squeeze y axis 
if numel(coIx)<=3
  niceyax(2);
else
  set(sp.selectCutout.axH,'ylim',[yOffs*1.02 r.aPeak(coIx(1))*.21]);
end
set(sp.selectCutout.axH,'xlim',wp.timeAx_pts([1 end]));
hold off

%   switch ap.nDecayComponent
%     case 1
%       f=@(t,a1,tau1) a1*exp(-t/tau1);
%     case 2
%       f=@(t,a1,tau1,a2,tau2) a1*exp(-t/tau1) + a2*exp(-t/tau2);
%   end


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
% §§ when introducing new parameters with time as unit make sure conversion
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
'10-90% rise time',...
'fit quality',...
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
  out=repmat(nan,numel(r.base),numel(item));
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
      case '10-90% rise time'
        out(:,g)=r.tRise'*tConvFac;
      case 'fit quality'
        out(:,g)=r.qFit';
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
