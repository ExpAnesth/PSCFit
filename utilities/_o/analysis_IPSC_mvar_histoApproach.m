function analysis_IPSC_mvar
% ** function analysis_IPSC_mvar
% multivariate analysis of IPSCs - histogram approach


% -------------------- PREPARATIONS ---------------------------------------
compName=lower(getenv('computername'));
switch compName
  case {'hh-i7'}
    dataPath='d:\_data\otc_ctx\ACh\AChBlockIPSC\';
    plotPath='d:\_data\otc_ctx\ACh\AChBlockIPSC\';
  case {'hh64','hh-i5'}
    dataPath='e:\_data\otc_ctx\ACh\AChBlockIPSC\';
    plotPath='e:\_data\otc_ctx\ACh\AChBlockIPSC\';
  case {'eval_lmb'}
    dataPath='h:\_data\otc_ctx\ACh\AChBlockIPSC\';
    plotPath='h:\_data\otc_ctx\ACh\AChBlockIPSC\';
  otherwise
    error('machine not defined');
end


% subdir
dSubDir='ACh_ACh_Block\Figs\';
% name of data file
dFn='AChAChBlock.mat';
% indexes into PSCRMN, corresponding to drug conditions listed in
% ds.indepPar and ds.indepParLabel
indepParIx=[1 2 3];
indepParTag={'ACh0','ACh+','ACh-'};


% subdir
dSubDir='Ctrl-ACh-Block\Figs\';
% name of data file
dFn='AChBlock.mat';
% indexes into PSCRMN, corresponding to drug conditions listed in
% ds.indepPar and ds.indepParLabel
indepParIx=[1 2 3];
indepParTag={'ACh0','ACh+','ACh-'};

% % subdir
% dSubDir='ACh-Dia-Block\Figs\';
% % name of data file
% dFn='AChDiaBlock.mat';
% % indexes into PSCRMN, corresponding to drug conditions listed in
% % ds.indepPar and ds.indepParLabel
% indepParIx=[1 2 3];
% indepParTag={'ACh0','ACh+','ACh-'};
 
% % subdir
% dSubDir='ACh-Zolpi200-Block\Figs\';
% % name of data file
% dFn='AChZolpi200Block.mat';
% % indexes into PSCRMN, corresponding to drug conditions listed in
% % ds.indepPar and ds.indepParLabel
% indepParIx=[1 2 3];
% indepParTag={'ACh0','ACh+','ACh-'};

% 
% % subdir
% dSubDir='ACh-Zolpi-Block\Figs\';
% % name of data file
% dFn='AChZolpiBlock.mat';
% % indexes into PSCRMN, corresponding to drug conditions listed in
% % ds.indepPar and ds.indepParLabel
% indepParIx=[1 2 3];


nIdepPar=numel(indepParIx);
% indexes into PSCRMN for statistical comparisons
idepStatsIx=[2 3];

% ** this variable determines which kind of IPSC parameters are analyzed:
% - 'detected': from all detected IPSCs, including those which could be fitted
% - 'fitted': from fitted IPSCs only
anParType='detected';

switch anParType
  case 'fitted'
    % - parameters to analyze
    % - description including units
    % - bins for 2D hist (at the same time, these define the limits of
    % axes in all plots)
    % - ticks for axes
    anPar={...
      'tDecay', '{\tau}_{decay} (ms)',2.^[1:.1:7],[10 100];
      'amp', 'peak amplitude (pA)',2.^[4:.1:11],[10 100 1000];
      %   'tRise20_80' , '20-80% rise time (ms)',2.^[-3.5:.1:3.4],[1 10] ;
      %    'width', 'width at half-ampl. (ms)', ;
      %  'chargePPsc', 'charge per PSC (pA*ms)',2.^[-2:.1:3],[.1 1]; ...
      };
  case 'detected'
    anPar={...
      'allTRise20_80' , '20-80% rise time (ms)',2.^[-3.5:.1:3.4],[1 10] ;
      'allAmp', 'peak amplitude (pA)',2.^[3.5:.1:11],[10 100 1000];
      };
end

nAnPar=size(anPar,1);

doPCA=false;
doNormalize=true;
doLogTransform=false;
% set to zero for no bootstrapping 
nBoot=00;
doHistStats=nBoot>0;
% if true, surface/contour plots will have same color scaling across drug
% conditions
doSameColorScale=true;
% if true, histogram data will be scaled according to the number of
% underlying events (wich is independent of parameter doSameColorScale)
doHistScale=true;

if doPCA
  % define number of PCs
  pcStr={'PC1','PC2','PC3'};
  nPC=numel(pcStr);
  nAnPar=nPC;
  % number of histogram bins in each dimension
  nHistBin=80;
end

% graphics
% define colors corresponding to ACh status 
fac=graphicsDef_ACh;
fIx=find(strcmp('ACh status',{fac.name}));
pCol=[];
for dcIx=1:nIdepPar
  pCol(dcIx,:)=fac(fIx).color(strcmp(indepParTag{dcIx},fac(fIx).levelName),:);
end

labelscale('fontSz',8,'scaleFac',1,'lineW',1,'markSz',3); 

% - which drug concs? More generally, which subset of data to use for
% defining PCs? Maybe data from all experiments, but these experiments
% plotted and analyzed separately?

% do it
load([dataPath dSubDir dFn]);
% index into 'slices' of PSCR (=parameters)
[~,~,parIx]=intersect(anPar(:,1)',depPar,'stable');
% number of nonredundant combinations of PCs/parameters to analyze and
% plot
[c,nc]=combin(nAnPar,'autoC',0);
% index into frequency parameters, for which there is only one value per
% cell, and which will be needed for display purposes
[~,~,freqParIx]=intersect({'freq','freqFit'},depPar,'stable');
% index into 'weight' parameter 
% ** make sure sensible things will be done if detected (vs fitted)
% parameters are loaded, for which this one is useless
switch anParType
  case 'detected'
    [~,~,weightParIx]=intersect({'allAmp'},depPar,'stable');
  case 'fitted'
    [~,~,weightParIx]=intersect({'chargePPsc'},depPar,'stable');
end

% 2D histogram templates
for pci=1:nc
  if doPCA
    histTemplate{pci}=nan(nHistBin,nHistBin,nBoot+1,nIdepPar);
    bootHistTemplate{pci}=nan(nHistBin,nHistBin);
  else
    histTemplate{pci}=nan(numel(anPar{c(pci,1),3}),numel(anPar{c(pci,2),3}),nBoot+1,nIdepPar);
    % template for bootstrapped data of compared drug conditions:
    bootHistTemplate{pci}=nan(numel(anPar{c(pci,1),3}),numel(anPar{c(pci,2),3}));
  end
end

% 2d Gaussian for smoothing
gauss2d=Gauss2D;

extractPar=nan(size(PSCR,1),1);
% loop over experiments
for g=1:size(PSCR,1)
  % -----------------------------------------------------------------------
  %                       COMPUTATIONS
  % -----------------------------------------------------------------------
  % collect data:
  dAll=[];
  groupAll=[];
  weightAll=[];
  nObsAll=nan(1,nIdepPar);
  % loop over concs
  for dcIx=1:nIdepPar
    % - observations in rows, parameters in columns
    curD=cat(2,PSCR{g,indepParIx(dcIx),parIx});
    curW=cat(2,PSCR{g,indepParIx(dcIx),weightParIx});
    % - get rid of entries with any nan or nonpositive number
    badIx=any(isnan(curD) | curD<=0,2);
    curD(badIx,:)=[];
    curW(badIx,:)=[];
    % - number of observations
    nObsAll(dcIx)=size(curD,1);
    % - group ID
    curGroup=repmat(indepParIx(dcIx),[nObsAll(dcIx) 1]);
    % - log transform?
    if doLogTransform
      curD=log(curD);
    end
    dAll=cat(1,dAll,curD);
    groupAll=cat(1,groupAll,curGroup);
    weightAll=cat(1,weightAll,curW);
  end
  
  
  % PCA:
  if doPCA
    pcaGroup=groupAll;
    [pc,td,v,ve,delIx]=PCexplore(dAll,'tag',groupAll,'nPC',nPC,'normalize',doNormalize,'plotType',[]);
    pcaGroup(delIx)=[];
    % now copy values to variables to be used further below
    d=td;
    group=pcaGroup;
    % define axis labels
    axLabel=pcStr;
  else
    d=dAll;
    group=groupAll;
    axLabel=anPar(1:nAnPar,2);
  end
  
  
  % 2D histograms & analysis:
  % each element of histPony contains the 2D histograms of the data
  % bins par1 | bins par2 | boot index | drug cond
  histPony=histTemplate;
  % bins par1 | bins par2 
  histHorse=bootHistTemplate;
  % - loop over parameter combinations (same number as PC)
  for pci=1:nc
    
    % - loop over concs
    for dcIx=1:nIdepPar
      curIx=group==indepParIx(dcIx);
      % ** retain only two columns as these are under scrutiny here
      curD=d(curIx,c(pci,:));
      nObs=nObsAll(dcIx);
      if doHistScale
        hScaleFac=1/nObs;
      else
        hScaleFac=1;
      end
      
      if doPCA
        % 2D histogram
        [tmpH,hCenter]=hist3(curD,nHistBin*[1 1])*hScaleFac;
        % use the occasion to fill the bin slots of anPar
        anPar(c(pci,:),3)=hCenter';
      else
        tmp=numel(find(curD(:,1)>anPar{c(pci,1),3}(end)))/nObs;
        if tmp>.01
          warning([num2str(tmp*100) ' % of events not covered by bin edges of par ' anPar{c(pci,1),1}]);
        end
        tmp=numel(find(curD(:,2)>anPar{c(pci,2),3}(end)))/nObs;
        if tmp>.01
          warning([num2str(tmp*100) ' % of events not covered by bin edges of par ' anPar{c(pci,2),1}]);
        end
        
        % § 2D histogram
        tmpH=hist3(curD,'edges',anPar(c(pci,:),3))*hScaleFac;
        % tmpH=weighHist2D(curD,weightAll,anPar{c(pci,1),3},anPar{c(pci,2),3})*hScaleFac;
        
      end
      
      % smooth by convolution with 2D Gaussian:
      tmpH=conv2(tmpH,gauss2d,'same');
      histPony{pci}(:,:,1,dcIx)=tmpH;
      
      if doHistStats
        % generate bootstrapped versions of data in loop, generate hist,
        % normalize to number of events, transform to row array and stuff
        % in histPony
        for bi=1:nBoot
          bootCurD=datasample(curD,nObs);
          if doPCA
            % 2D histogram
            tmpH=hist3(bootCurD,hCenter)*hScaleFac;
          else
            tmpH=hist3(bootCurD,'edges',anPar(c(pci,:),3))*hScaleFac;
          end
          % smooth
          tmpH=conv2(tmpH,gauss2d,'same');
          histPony{pci}(:,:,bi+1,dcIx)=tmpH;
        end
      end
    end
    
    if doHistStats
      % pick pair to be compared statistically
      [~,tmpIx]=intersect(indepParIx,idepStatsIx,'stable');
      % compare UNSCALED histograms bin by bin
      x=histPony{pci}(:,:,2:end,tmpIx(1))*nObsAll(tmpIx(1));
      x=reshape(permute(x,[3 2 1]),[nBoot,numel(anPar{c(pci,1),3}) * numel(anPar{c(pci,2),3})]);
      
      y=histPony{pci}(:,:,2:end,tmpIx(2))*nObsAll(tmpIx(2));
      y=reshape(permute(y,[3 2 1]),[nBoot,numel(anPar{c(pci,1),3}) * numel(anPar{c(pci,2),3})]);
      
      % now compare bin by bin
      stats=mes(x,y,{'auroc','hedgesg'});
      % reshape and save stats values
      histHorse{pci}=reshape(stats.hedgesg,[numel(anPar{c(pci,2),3}),numel(anPar{c(pci,1),3})]);
    end
  end
  % -----------------------------------------------------------------------
  %                       PLOTTING
  % -----------------------------------------------------------------------
  
  fh=mkfig(g,'b');
  set(fh,'name',LISTEXP{g});
  colormap(coma('CubicL'));
  clf
  % ensure that plots are not too elongated
  nRow=max(3,nc);
  nCol=nIdepPar+3;

  mima=cellfun(@(x) x([1 end])',anPar(:,3),'UniformOutput',false);
  mima=cat(2,mima{:});

  % - loop over parameter combinations (same number as PC)
  for pci=1:nc
    % --- SCATTER PLOTS ---
    subplot(nRow,nCol,(pci-1)*nCol+1)
    hold on
    % raw points (plot for all drug conditions so that averages won't be overplotted)
    for dcIx=1:nIdepPar
      curIx=group==indepParIx(dcIx);
      ph=plot(d(curIx,c(pci,1)),d(curIx,c(pci,2)),'.');
      set(ph,'color',pCol(dcIx,:));
    end
    % - loop over drug conditions: averages & errors
    for dcIx=1:nIdepPar
      curIx=group==indepParIx(dcIx);
      if 1
        prc=prctile(d(curIx,:),[5 50 95]);
        av=prc(2,:);
        er=nan*prc([1 3],:);
      elseif 0
        av=nanmean(d(curIx,:));
        st=nanstd(d(curIx,:));;
        er=[av-st; av+st];
      else
        av=geomean(d(curIx,:));
        er=nan*[av; av];
      end
      
      lh=errorcross2(...
        [av(c(pci,1)),av(c(pci,2))],...
        [av(c(pci,1))-er(1,c(pci,1))  av(c(pci,2))-er(1,c(pci,2))],...
        [er(2,c(pci,1))-av(c(pci,1))  er(2,c(pci,2))-av(c(pci,2))],...
        'color',pCol(dcIx,:)*.7,'linewidth',2.0);
      
      ph=plot(av(c(pci,1)),av(c(pci,2)),'o');
      set(ph,'color',pCol(dcIx,:)*.7,'markerfacecolor',pCol(dcIx,:),'markersize',8);
      
      if dcIx==nIdepPar
        if doPCA
          axis tight
        else
          set(gca,'xscale','log','yscale','log');
          axis([mima(:,c(pci,1))' mima(:,c(pci,2))'])
        end
        grid on
        xlabel(axLabel{c(pci,1)});
        ylabel(axLabel{c(pci,2)});
      end
    end
    
    % --- SURFACE PLOTS ---
    if doSameColorScale
      tmp=histPony{pci}(:,:,1,:);
      clim=prctile(tmp(:),[.5 99.9]);
    end
    for dcIx=1:nIdepPar
      sph=subplot(nRow,nCol,(pci-1)*nCol+dcIx+1);
      if any(any(isfinite(histPony{pci}(:,:,1,dcIx))))
        % unfortunately, grid does not show up in surface plot viewed from
        % above
        % h=surf(anPar{c(pci,1),3},anPar{c(pci,2),3},histPony{pci}(:,:,1,dcIx)');
        % view(0,90);
        [~,h]=contourf(anPar{c(pci,1),3},anPar{c(pci,2),3},histPony{pci}(:,:,1,dcIx)',40);
        set(h,'linestyle','none');
        if doSameColorScale
          set(gca,'cLim',clim);
        end
        if ~doPCA
          set(gca,'xscale','log','yscale','log');
          set(gca,'xtick',anPar{c(pci,1),4},'ytick',anPar{c(pci,2),4})
        end
        grid on
        grid minor
        sph.GridAlpha=.8;
        sph.GridLineStyle='--';
        % same axis limits as corresponding scatter plot
        axis([mima(:,c(pci,1))' mima(:,c(pci,2))'])
        % info on IPSC freq
        smarttext(['f=' num2str(PSCR{g,indepParIx(dcIx),freqParIx(1)},2) '/'...
          num2str(PSCR{g,indepParIx(dcIx),freqParIx(2)},2) ' Hz'],.0,.8,'color','w','fontweight','bold');
        xlabel(axLabel{c(pci,1)});
      end
      title(indepParTag{dcIx},'fontsize',12);
    end
    
  end
  
  
  % difference of UNSCALED histograms
  sph=subplot(nRow,nCol,(pci-1)*nCol+nIdepPar+2);
  colormap(sph,coma('bluered'));
  [~,tmpIx]=intersect(indepParIx,idepStatsIx,'stable');
  tmpH=histPony{pci}(:,:,1,tmpIx(1))*nObsAll(tmpIx(1)) - ...
    histPony{pci}(:,:,1,tmpIx(2))*nObsAll(tmpIx(2));
  cc=corr2(histPony{pci}(:,:,1,tmpIx(1)),histPony{pci}(:,:,1,tmpIx(2)))
  %   % Bhattacharyya distance:
  %   bd=sum(sum(sqrt(histPony{pci}(:,:,1,tmpIx(1)).*histPony{pci}(:,:,1,tmpIx(2)))))
  
  if any(any(isfinite(tmpH)))
    % symmetric color limits
    cl=max(abs(tmpH(:)))*[-1 1];
    [~,h]=contourf(anPar{c(pci,1),3},anPar{c(pci,2),3},tmpH',40);
    set(h,'linestyle','none');
    set(gca,'clim',cl);
    if ~doPCA
      set(gca,'xscale','log','yscale','log');
      set(gca,'xtick',anPar{c(pci,1),4},'ytick',anPar{c(pci,2),4})
    end
    grid on
    grid minor
    sph.GridAlpha=.8;
    sph.GridLineStyle='--';
    % same axis limits as corresponding scatter plot
    axis([mima(:,c(pci,1))' mima(:,c(pci,2))'])
    smarttext(['cc=' num2str(cc)],.0,.8,'color','k','fontweight','bold');
    colorbar('location','east')
  end

  
  extractPar(g)=cc
  

  
%   % stats
%   sph=subplot(nRow,nCol,(pci-1)*nCol+nIdepPar+3);
%   colormap(sph,coma('bluered'));
%   tmpH=histHorse{pci}';
%   if any(any(isfinite(tmpH)))
% %     % color limits
% %     cl=[0 1];
%     % symmetric color limits
%     cl=max(abs(tmpH(:)))*[-1 1];
%     
%     [~,h]=contourf(anPar{c(pci,1),3},anPar{c(pci,2),3},tmpH',40);
%     set(h,'linestyle','none');
%     set(gca,'clim',cl);
%     if ~doPCA
%       set(gca,'xscale','log','yscale','log');
%       set(gca,'xtick',anPar{c(pci,1),4},'ytick',anPar{c(pci,2),4})
%     end
%     grid on
%     grid minor
%     sph.GridAlpha=.8;
%     sph.GridLineStyle='--';
%     % same axis limits as corresponding scatter plot
%     axis([mima(:,c(pci,1))' mima(:,c(pci,2))'])
%     colorbar
%   end
    
  
  
%   if doPCA
%     subplot(2,4,4)
%     % plot PC
%     ph=plot(pc,'-o');
%     nicexyax;
%     ylabel('weight')
%     legend(ph,pcStr)
%     grid on
%     set(gca,'xtick',1:nAnPar,'xticklabel',anPar(:,1))
%     rotateXLabels(gca,60)
%     xlabel('parameter');
%   end

  drawnow
  pause(.05)
  
end
  
  
function g=Gauss2D
% support up to 3 sigma, sigma corresponding to 2 bins
mu = [0 0];
Sigma = [1 0;0 1];
% step size 1: sigma=1
x1 = -3:.75:3; x2 = -3:.75:3;
[X1,X2] = meshgrid(x1,x2);
g = mvnpdf([X1(:) X2(:)],mu,Sigma);
% normalize so that elements sum to 1
g=g/sum(g);
g = reshape(g,length(x2),length(x1));


function [h,varargout]=weighHist2D(d,weight,bin1,bin2)
% produce 2D histogram of weighted data (h) and, optionally, a cell array
% of indexes to corresponding elements

n1=numel(bin1);
n2=numel(bin2);
bin1(end+1)=inf;
bin2(end+1)=inf;
h=zeros(n1,n2);
if nargout>1
  doCollectIdx=true;
  ih=cell(n1,n2);
else
  doCollectIdx=false;
end

for b1Ix=1:n1
  for b2Ix=1:n2
    ix=d(:,1)>=bin1(b1Ix) & d(:,1)<bin1(b1Ix+1) & d(:,2)>=bin2(b2Ix) & d(:,2)<bin2(b2Ix+1);
    h(b1Ix,b2Ix)=sum(weight(ix));
    if doCollectIdx
      ih{b1Ix,b2Ix}=find(ix);
    end
  end
end
    