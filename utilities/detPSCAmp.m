function [amp,tRise,varargout]=detPSCAmp(dExc,ddExc,tZeroIx,si,tsl,varargin)
% ** function [amp,tRise,varargout]=detPSCAmp(dExc,ddExc,tZeroIx,thresh,si,tsl,varargin)
% estimates amplitudes and rise times of phasic postsynaptic currents
% (PSCs) from/of their rise phase. The input arguments to detPSCAmp must
% have been produced in the following way:
% 1. Lowpass filter the raw data trace, eliminating any phase lag (e.g.
% using filtfilt)
% 2. Run these data through a differentiator filter, e.g. via
% designfilt('differentiatorfir',...). NOTES: a) make sure that the corner
% frequency of the differentiator lowpass is way higher than the corner
% frequency of the lowpass filter in 1.; b) eliminate phase delays imposed
% by the differentiator filter (which can be found via grpdelay). Function
% differfi accomplishes most of this.
% 3. In the lowpass-filtered derivative of the data produced in 2., detect
% rising phases of PSCs via a fixed threshold (which must also be furnished
% as input argument thresh)
% 4. Based on the time stamps obtained in 3. (input argument tsl) produce
% cutouts of PSCs from both original trace and derivative, (input variables
% dExc and ddExc, respectively). Reasonable cutout windows would be [-5
% 10] ms, depending of course on the kinetics of the rise time. Input
% argument tZeroIx is the index to the data points in the excerpt
% respresenting time zero.
% The code assumes positive-going events and will automatically invert the
% data if they are negative-going. It has been designed for bursts of EPSCs
% as typically occurring during full activity in interconnected networks,
% but should of course also work on 'classical' current traces with PSCs
% occurring mostly in isolation.
% Results of the computations can be visualized in two different ways. If
% optional input argument 'plotSingle' is true, detailed results of the
% computations will be shown for each cutout. If input argument 'nPlotEv'
% is a positive integer a plot of parts of the raw data and the computed
% amplitudes will be produced.
% Optional input arguments (see below) must be specified as parameter/value
% pairs, e.g. as in
%          detPSCAmp(...,'nPlotEv',100);
%
%                         >>> INPUT VARIABLES >>>
%
% NAME           TYPE/DEFAULT       DESCRIPTION
% dExc           2D array           PSC excerpts (cutouts) from raw data
% ddExc          2D array           PSC excerpts from differentiated data 
% tZeroIx        scalar             index to first post-event point in cutouts
% si             scalar             sampling interval of data in µs
% tsl            1D array           time stamp list of events
% -------------------------------------------------------------------------
% specify following optional input args for overview plots only
% -------------------------------------------------------------------------
% fh             figure handle      handle to figure to be used for
%                                    displaying results (a new figure will 
%                                    be created if not specified)
% d              vector             raw time series from which dExc was 
%                                    produced
% thresh         scalar             threshold that was used for detection
%                                    of PSCs
% nPlotEv        scalar, 0          number of PSCs to plot in overview 
%                                    figure (set to zero or [] for no plot)
% plotSingle     logical, false     if true, detailed results will be shown
%                                    interactively for each cutout (halts 
%                                    computations!)
%
%                         <<< OUTPUT VARIABLES <<<
% NAME           TYPE/DEFAULT           DESCRIPTION
% amp            vector                 estimated PSC amplitudes (original
%                                        unit)
% tRise          vector                 estimated PSC rise time (ms)
% varargout{1}   vector                 estimated PSC base amplitudes
% varargout{2}   vector                 estimated PSC peak amplitudes


% TO DO:

% ----- default values & varargin -----
fh=[];
d=[];
thresh=nan;
nPlotEv=0;
plotSingle=false;

pvpmod(varargin);
disp(['** ' mfilename])

% -----------------------------------------------------------------------------
%                          I. PRELIMINARIES
% -----------------------------------------------------------------------------

% size of things
[n1,n2]=size(dExc);
nTs=size(tsl,1);
% some plausibility checks of input args
if tZeroIx<=0
  error('input arg tZero must be strictly positive');
elseif tZeroIx>=n1
  error('input arg tZero must be smaller than the length of the cutouts');
end
if ~isequal([n1,n2],size(ddExc))
  error('input args dExc and ddExc must be of same size')
end
if n2~=nTs
  error('number of cutouts and number of time stamps do not match')
end

% check polarity of data
polarity=sign(diff(diff(prctile(ddExc(:),[1 50 99]))));
if ~polarity
  error('funny data')
elseif polarity<0
  % invert
  disp('inverting data');
  dExc=dExc*polarity;
  ddExc=ddExc*polarity;
  thresh=thresh*polarity;
  if ~isempty(d)
    d=d*polarity;
  end
end

% check if plots are to be produced and all necessary input was provided
if ~isempty(nPlotEv) && nPlotEv>0
  if ~isempty(d) && ~isempty(si)
    doPlot=true;
    % place a few variables in struct gr to be used with guidata
    gr.nPlotEv=nPlotEv;
    gr.thresh=thresh;
    % x coordinate of event clicked upon in overview plot - set to empty
    % for random selection
    gr.x=[];
    % create figure and put out info so the user knows it's working
    if isempty(fh)
      fh=mkfig([],'b');
      set(fh,'Tag','evAmpFigure');
    else
      figure(fh)
    end
    % some default graphics settings:
    % - background colors of axes for better visibility
    gr.bgCol=[.9 .9 .9];
    gr.mCol1=[0 .5 0];
    gr.mCol2=[.8 .1 .6];
    gr.mSz=2;
    
    clf
    gr.sphOv=subplot(2,1,1);
    % § look deeper into graphics to set axes positions better
    rexy('ax',gca,'xfac',1.2,'yfac',1.05);
    gr.th=smarttext(['computing amplitudes of ' int2str(n2) ' events...'],.1,.5,'Fontsize',20);
    guidata(fh,gr);
    drawnow
  else
    doPlot=false;
    warning('input variables d and si are needed for results plots');
  end
else
  doPlot=false;
end

% -----------------------------------------------------------------------------
%                          II. COMPUTATIONS
% -----------------------------------------------------------------------------
% we need to compute all peaks in the partly preconditioned raw data
peakDExc=evdeal(dExc,'idx',{'allpeaks'});
% same for derivative
peakDdExc=evdeal(ddExc,'idx',{'allpeaks'});
% base of each detected event: going BACKWARDS from tZero, the event
% detection time stamp, it is either of
% - the first trough of dExc (dExc=cutouts of data)
% - the first neg->pos zero line crossing of ddExc 
%   (ddExc=cutouts of differentiated data)
% - the first trough of ddExc, whichever comes first.
% So, find them:
% - in dExc pick last trough before and including tZero
lastNegPeakIx=cellfun(@helperfun1,repmat({tZeroIx},1,n2),peakDExc.negPeakT);
lastNegPeakIx=lastNegPeakIx';
% - in ddExc up to ts, find indexes to negative points
[rw,cl]=find(ddExc(1:tZeroIx,:)<=0);
% - warn if there are cutouts for which none exist and set the
%   corresponding index to nan
[rw,cl]=rectifyCl(rw,cl,n2,{'pre','left'});
% - for each cutout, pick the last of the found pre-tZero indexes
baseZeroCrossIx=accumarray(cl,rw,[],@max);
% - in ddExc pick last trough before tZero
lastNegPeakDerivIx=cellfun(@helperfun1,repmat({tZeroIx},1,n2),peakDdExc.negPeakT);
lastNegPeakDerivIx=lastNegPeakDerivIx';
% - finally, of the three indexes found above, pick the one closest to tZero
baseIx=max([baseZeroCrossIx lastNegPeakIx lastNegPeakDerivIx],[],2);

% peak of each detected event: going FORWARDS from tZero, the event
% detection time stamp, it is either of
% - the first peak of dExc
% - the first pos->neg zero line crossing of ddExc
% - the first trough of ddExc, whichever comes first. So, find them:
% - in dExc starting with ts, find first peak
firstPosPeakIx=cellfun(@helperfun2,repmat({tZeroIx},1,n2),peakDExc.posPeakT);
firstPosPeakIx=firstPosPeakIx';
% - in ddExc starting from ts, find indexes to negative points
[rw,cl]=find(ddExc(tZeroIx:end,:)<0);
% - warn if there are cutouts for which none exist and set the
%   corresponding index to nan
[rw,cl]=rectifyCl(rw,cl,n2,{'post','right'});
% - for each cutout, pick the last point before first negative point
peakZeroCrossIx=tZeroIx + accumarray(cl,rw,[],@min) - 2;
% - in ddExc starting with ts, find first trough, which may be the only
% conceivable peak marker for a PSC with an immediate follower
firstNegPeakDerivIx=cellfun(@helperfun2,repmat({tZeroIx},1,n2),peakDdExc.negPeakT);
firstNegPeakDerivIx=firstNegPeakDerivIx';
% - finally, of the indexes found above, pick the one closest to tZero
peakIx=min([peakZeroCrossIx firstPosPeakIx firstNegPeakDerivIx],[],2);

diffBaseIx=diff(baseIx+tsl);
diffPeakIx=diff(peakIx+tsl);
% - events with identical peak and baseline indices to be eliminated: leave
% first)
identPeakAndBaseIx=[false; ~diffPeakIx & ~diffBaseIx];
% - events with identical peak but different base line indices to be eliminated
identPeakIx=[false; ~diffPeakIx & diffBaseIx];
% - events with identical base line but different peak indices to be eliminated
identBaseIx=[false; diffPeakIx & ~diffBaseIx];

if sum(identPeakAndBaseIx)
  disp(['eliminating ' int2str(sum(identPeakAndBaseIx)) ' events with overlapping base and peak indices'])
end
if sum(identPeakIx)
  disp(['eliminating ' int2str(sum(identPeakIx)) ' events with overlapping peak indices'])
end
if sum(identBaseIx)
  disp(['eliminating ' int2str(sum(identBaseIx)) ' events with overlapping base indices'])
end
delIx=identPeakAndBaseIx | identPeakIx | identBaseIx;

% use these indexes to compute base line and peak...
baseVal=dExc(sub2ind([n1,n2],baseIx,(1:n2)'));
peakVal=dExc(sub2ind([n1,n2],peakIx,(1:n2)'));
% but then set indexes and values of all events deemed bad to nan
baseIx(delIx)=nan;
baseVal(delIx)=nan;
peakIx(delIx)=nan;
peakVal(delIx)=nan;

% amplitudes 
amp=abs(peakVal-baseVal);
% rise time: 
% - upsample to 100 kHz, subtracting baseline values beforehand to minimize
% edge transients
uSi=10;
uDExc=resample(dExc-baseVal',si,uSi);
uBaseIx=(baseIx-1)*(si/uSi)+1;
uPeakIx=(peakIx-1)*(si/uSi)+1;
% - as base 'line' and peak are not restrained within specific windows we
% have to loop over excerpts
tRise=nan(n2,2);
badCount=0;
for k=setdiff(1:n2,find(delIx))
  curExc=uDExc(uBaseIx(k):uPeakIx(k),k);
  tmpTr=(find(curExc>=amp(k)*.9,1) - find(curExc>=amp(k)*.1,1))*(uSi/1000);
  tmpTr2=(find(curExc>=amp(k)*.8,1) - find(curExc>=amp(k)*.2,1))*(uSi/1000);
  if ~isempty(tmpTr)
    tRise(k,1)=tmpTr;
  else
    badCount=badCount+1;
  end
  if ~isempty(tmpTr2)
    tRise(k,2)=tmpTr2;
  end
end

if badCount>0
  warning(['rise time of ' int2str(badCount) ' events could not be determined'])
end
if nargout>=3
  varargout{1}=baseVal;
end
if nargout>=4
  varargout{2}=peakVal;
end


% -----------------------------------------------------------------------------
%                          III. PLOTS
% -----------------------------------------------------------------------------
if doPlot
  plotTsl=discrete2cont(tsl,si/1000);
  % index to subset of events to be plotted
  tmpPlotTsIx=1:min(nPlotEv,n2);
  % focus on figure
  figure(fh);
  % remove info text
  delete(gr.th);
  % time in real units
  tmpT=discrete2cont(tsl(tmpPlotTsIx(1)):tsl(tmpPlotTsIx(end)),si/1000);
  % plot raw data and lines depicting computed event amplitudes of up to
  % nPlotEv events, starting from beginning of time series
  plot(tmpT,d(tsl(tmpPlotTsIx(1)):tsl(tmpPlotTsIx(end))));
  niceyax;
  hold on
  lh=line(cat(1,(plotTsl(tmpPlotTsIx))',...
    (plotTsl(tmpPlotTsIx))'+discrete2cont(peakIx(tmpPlotTsIx)'-tZeroIx,si/1000)),...
    [baseVal(tmpPlotTsIx) peakVal(tmpPlotTsIx)]','color','m','linewidth',1.5);
  set(gr.sphOv,'color',gr.bgCol);
  xlabel('time (ms)')
  % callback for lines marking rise phases
  set(lh,'ButtonDownFcn',{@plotSelectedExcerpts,fh,plotTsl,tZeroIx,dExc,ddExc,baseIx,baseVal,peakIx,peakVal});
  % plot a handful of excerpts 
  plotSelectedExcerpts([],[],fh,plotTsl,tZeroIx,dExc,ddExc,baseIx,baseVal,peakIx,peakVal);
end
drawnow
if plotSingle
  figure
  pscInspector(dExc,ddExc,1:size(dExc,2),lastNegPeakIx,firstPosPeakIx,baseIx,peakIx,baseZeroCrossIx,lastNegPeakDerivIx,peakZeroCrossIx,firstNegPeakDerivIx)
end

% ========================== LOCAL FUNCTIONS ==============================
function plotSelectedExcerpts(src,~,fh,plotTsl,tZeroIx,dExc,ddExc,baseIx,baseVal,peakIx,peakVal)
% retrieve data
gr=guidata(fh);
if ~isempty(src)
  % set current x coordinate
  gr.x=get(src,'x');
  gr.x=gr.x(1);
  % store data
  guidata(fh,gr);
end

n1=size(dExc,1);
nPlotCol=5;

if isempty(gr.x)
  plotEvIx=randperm(gr.nPlotEv,nPlotCol);
else
  [~,plotEvIx]=min(abs(gr.x-plotTsl));
  plotEvIx=(-floor(nPlotCol/2):floor(nPlotCol/2))+plotEvIx-1;
  plotEvIx=plotEvIx+max(1,-plotEvIx(1));
  %§
end

for spIx=1:numel(plotEvIx)
  evIx=plotEvIx(spIx);
  sph=subplot(2,nPlotCol,nPlotCol+spIx);
  cla
  % offset for raw data traces: mean
  tmpOffs=mean(dExc(:,evIx));
  % plot cutouts from raw data WITH offset
  ph=plot(dExc(:,evIx)-tmpOffs,'-o');
  set(ph,'markersize',gr.mSz);
  hold on
  % plot cutouts from differentiated data WITHOUT offset
  ph=plot(ddExc(:,evIx),'-o');
  set(ph,'markersize',gr.mSz);
  % mark base and peak
  ph=plot(baseIx(evIx),baseVal(evIx)-tmpOffs,'o');
  set(ph,'markersize',gr.mSz*3,'linewidth',1,'color',gr.mCol1);
  ph=plot(peakIx(evIx),peakVal(evIx)-tmpOffs,'mo');
  set(ph,'markersize',gr.mSz*3,'linewidth',1,'color',gr.mCol2);
  niceyax;
  smarttext(['#' int2str(evIx)],.1,.85);
  yl=get(gca,'ylim');
  % line indicating detection tZero
  line((tZeroIx)*[1 1],yl,'color','k','linestyle','--');
  % line for threshold
  line([1 n1],gr.thresh*[1 1],'color','k','linestyle',':');
  % zero line 
  line([1 n1],[0 0],'color','k','linestyle','--');
  % bg col
  set(sph,'color',gr.bgCol);
  set(gca,'xtick',[]);
end

function outTIx=helperfun1(tZeroIx,tIx)
ix=find(tIx<=tZeroIx,1,'last');
if ~isempty(ix)
  outTIx=tIx(ix);
else
  outTIx=nan;
end

function outTIx=helperfun2(tZeroIx,tIx)
ix=find(tIx>=tZeroIx,1,'first');
if ~isempty(ix)
  outTIx=tIx(ix);
else
  outTIx=nan;
end

function [rw,cl]=rectifyCl(rw,cl,n2,strCell)
clApp=setdiff((1:n2)',cl);
if ~isempty(clApp)
  warning(['in ' int2str(numel(clApp)) ' cutout(s) of differentiated data no ' strCell{1} '-detection zero-line crossing exists - consider extending ' strCell{2} ' border of cutout interval']);
  cl=cat(1,cl,clApp);
  rw=cat(1,rw,clApp*nan);
end

function pscInspector(dExc,ddExc,evIx,lastNegPeakIx,firstPosPeakIx,baseIx,peakIx,baseZeroCrossIx,lastNegPeakDerivIx,peakZeroCrossIx,firstNegPeakDerivIx)
% loop over all excerpts to see in detail how algorithm picks base and peak
for ix=evIx
  subplot(1,2,mod(ix-1,2)+1)
  cla
  yyaxis left
  plot(dExc(:,ix));
  hold on;
  plot(lastNegPeakIx(ix),dExc(lastNegPeakIx(ix),ix),'go');
  plot(firstPosPeakIx(ix),dExc(firstPosPeakIx(ix),ix),'mo');
  if isfinite(baseIx(ix))
    plot(baseIx(ix),dExc(baseIx(ix),ix),'ko','markerfacecolor','g');
  end
  if isfinite(peakIx(ix))
    plot(peakIx(ix),dExc(peakIx(ix),ix),'ko','markerfacecolor','m');
  end
  set(gca,'xgrid','on')
  niceyax
  hold off
  
  yyaxis right
  plot(ddExc(:,ix));
  hold on;
  plot(baseZeroCrossIx(ix),ddExc(baseZeroCrossIx(ix),ix),'go');
  if isfinite(lastNegPeakDerivIx(ix))
    plot(lastNegPeakDerivIx(ix),ddExc(lastNegPeakDerivIx(ix),ix),'go');
  end
  plot(peakZeroCrossIx(ix),ddExc(peakZeroCrossIx(ix),ix),'mo');
  if isfinite(firstNegPeakDerivIx(ix))
    plot(firstNegPeakDerivIx(ix),ddExc(firstNegPeakDerivIx(ix),ix),'mo');
  end
  
  niceyax
  line([1 size(dExc,1)],[0 0],'color','k','linestyle',':')
  hold off 
  title(['event # ' int2str(ix) '; key to proceed'])
  pause
end