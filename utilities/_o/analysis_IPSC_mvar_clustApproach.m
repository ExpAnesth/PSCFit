function analysis_IPSC_mvar
% ** function analysis_IPSC_mvar
% multivariate analysis of IPSCs - clustering approach

% to do:
% MAJOR
% - different clustering approach: cluster each drug condition separately,
% then merge clusters and count and compare number of events in the merged
% set of clusters as before?
% - use histogram/bin counts and Cramer's V for additional stats
% - to 
% SMALL FRY
% - what to do if a cluster is empty?
% - does bootstrapping make sense in its current form?


% -------------------- PREPARATIONS ---------------------------------------
compName=lower(getenv('computername'));
switch compName
  case {'hh-i7'}
    dataPath='d:\_data\otc_ctx\ACh\AChBlockIPSC\';
    plotPath='d:\_data\otc_ctx\ACh\AChBlockIPSC\';
  case {'hh64','hh-i5'}
    dataPath='e:\_data\otc_ctx\ACh\AChBlockIPSC\';
    plotPath='e:\_data\otc_ctx\ACh\AChBlockIPSC\';
  case {'eval-lmb'}
    dataPath='h:\_data\otc_ctx\ACh\AChBlockIPSC\';
    plotPath='h:\_data\otc_ctx\ACh\AChBlockIPSC\';
  case {'ll'}
    dataPath='d:\ll\Data\VoltageClamp\';
    plotPath='d:\ll\Data\VoltageClamp\';
  otherwise
    error('machine not defined');
end


% --- key user options
printas='-djpeg95';
printas=[];
% scaling factor for bubble plots: an amplitude of 1000 pA shall correspond
% to area 100;
ampScalFac=1/10;

% - clustering parameters:
% shall data be log-transformed for clustering? 
doLogTransform=true;
% if different from Nan, this fixes the number of clusters
nFixClust=5;
nMaxClust=8;
clustCol=[.4 .4 .4; 0 0 .5; .5 .3 0];
clustMeth='kmeans';
% clustMeth='hGauss';
% distance method for k-means
distType='sqeuclidean';
nReplic=10;
cluOpts=statset('display','final','UseParallel',false,'MaxIter',500);

% set to zero for no bootstrapping 
nBoot=0;
doBootStats=nBoot>0;
% if true, surface/contour plots will have same color scaling across drug
% conditions
doSameColorScale=true;
% if true, histogram data will be scaled according to the number of
% underlying events (wich is independent of parameter doSameColorScale)
doHistScale=true;

% indexes into PSCRMN to drug conditions to be used for statistical comparisons
drugStatsIx=[2 3];

% indexes into PSCRMN to drug condition(s) to be used for cluster generation
drugClusterIx=[2];

% subdir
dSubDir='ACh_ACh_Block\Figs\';
% name of data file
dFn='AChAChBlock.mat';
% indexes into PSCRMN, corresponding to drug conditions listed in
% ds.indepPar and ds.indepParLabel
drugIx=[1 2 3];
drugTag={'ACh0','ACh+','ACh-'};


% subdir
dSubDir='Ctrl-ACh-Block\Figs\';
% name of data file
dFn='AChBlock.mat';
% indexes into PSCRMN, corresponding to drug conditions listed in
% ds.indepPar and ds.indepParLabel
drugIx=[1 2 3];
drugTag={'ACh0','ACh+','ACh-'};

% % subdir
% dSubDir='ACh-Dia-Block\Figs\';
% % name of data file
% dFn='AChDiaBlock.mat';
% % indexes into PSCRMN, corresponding to drug conditions listed in
% % ds.indepPar and ds.indepParLabel
% drugIx=[1 2 3];
% drugTag={'ACh0','ACh+','ACh-'};
%  
% % subdir
% dSubDir='ACh-Zolpi200-Block\Figs\';
% % name of data file
% dFn='AChZolpi200Block.mat';
% % indexes into PSCRMN, corresponding to drug conditions listed in
% % ds.indepPar and ds.indepParLabel
% drugIx=[1 2 3];
% drugTag={'ACh0','ACh+','ACh-'};
% 
% % subdir
% dSubDir='ACh_TTX\Figs\';
% % name of data file
% dFn='AChTTX.mat';
% % indexes into PSCRMN, corresponding to drug conditions listed in
% % ds.indepPar and ds.indepParLabel
% drugIx=[1 2 3];
% drugTag={'ACh0','ACh+','ACh-'};


% % subdir
% dSubDir='ACh-Zolpi-Block\Figs\';
% % name of data file
% dFn='AChZolpiBlock.mat';
% % indexes into PSCRMN, corresponding to drug conditions listed in
% % ds.indepPar and ds.indepParLabel
% drugIx=[1 2 3];


% - parameters to load, check and possibly analyze
% - description including units
% - bins for 2D hist (at the same time, these define the limits of
%  acceptable values)
% - ticks for axes
% - type of par: - 'detected': from all detected IPSCs, including those which could be fitted
%                - 'fitted': from fitted IPSCs only

fullPSCPar={...
  'allTRise20_80' , '20-80% rise time (ms)',2.^[log2(.05):.11:2.6],[.1 1], 'detected';
  'allAmp', 'peak amplitude (pA)',2.^[log2(10):.13:11],[10 100 1000],      'detected';
  'tDecay', '{\tau}_{decay} (ms)',2.^[1:.1:7],[10 100],                    'fitted';
  'amp', 'peak amplitude (pA)',2.^[log2(10):.13:11],[10 100 1000],         'fitted';
  'tRise20_80' , '20-80% rise time (ms)',2.^[log2(.05):.11:2.6],[.1 1],    'fitted' ;
%  'chargePPsc', 'charge per PSC (pA*ms)',2.^[-2:.1:3],[.1 1], 'fitted'; ... % §§ check limits!
};

% *** parameters to be used for scatter and histogram plots and clustering
plotPSCParInd=[1 2];
plotPSCPar=fullPSCPar(plotPSCParInd,:);

labelscale('fontSz',8,'scaleFac',1,'lineW',1,'markSz',4); 

% -------------------------------------------------------------------------
% --------------------- run, rabbit, run! ---------------------------------
% -------------------------------------------------------------------------

% precompute some frequently needed variables
nDrugCondit=numel(drugIx);
nFullPSCPar=size(fullPSCPar,1);
nPlotPSCPar=size(plotPSCPar,1);
% number of nonredundant combinations of PCs/parameters to analyze and
% plot
[c,nc]=combin(nPlotPSCPar,'autoC',0);
% 2D histogram template
histTemplate=nan(numel(plotPSCPar{1,3}),numel(plotPSCPar{2,3}),nDrugCondit);
% 2d Gaussian for smoothing
gauss2d=Gauss2D;

% plots:
% define colors corresponding to ACh status 
fac=graphicsDef_ACh;
fIx=find(strcmp('ACh status',{fac.name}));
pCol=[];
for dcIx=1:nDrugCondit
  pCol(dcIx,:)=fac(fIx).color(strcmp(drugTag{dcIx},fac(fIx).levelName),:);
end
axLabel=plotPSCPar(1:nPlotPSCPar,2);

% LOAD DATA
load([dataPath dSubDir dFn]);

% indexes into data:
% - plot parameters (=slices of PSCR)
[~,~,parIx]=intersect(plotPSCPar(:,1)',depPar,'stable');
% - all parameters
[~,~,allTRise20_80Ix]=intersect({'allTRise20_80'},depPar,'stable');
[~,~,allAmpIx]=intersect({'allAmp'},depPar,'stable');
[~,~,tDecayIx]=intersect({'tDecay'},depPar,'stable');
[~,~,ampIx]=intersect({'amp'},depPar,'stable');
[~,~,tRise20_80Ix]=intersect({'tRise20_80'},depPar,'stable');

% - frequency parameters, for which there is only one value per cell, and
% which will be needed for display purposes
[~,~,freqParIx]=intersect({'freq','freqFit'},depPar,'stable');
% - time stamp lists:
[~,~,allTslIx]=intersect({'allTsl'},depPar,'stable');
[~,~,tslIx]=intersect({'tsl'},depPar,'stable');
% finally, the kinds of parameters:
detParInd=find(strcmp(fullPSCPar(:,5),'detected'));
fitParInd=find(strcmp(fullPSCPar(:,5),'fitted'));
[~,~,detParIx]=intersect(fullPSCPar(detParInd,1)',depPar,'stable');
[~,~,fitParIx]=intersect(fullPSCPar(fitParInd,1)',depPar,'stable');


% loop over experiments
for g=1:size(PSCR,1)
  % -----------------------------------------------------------------------
  % -----------------------------------------------------------------------
  %                       COMPUTATIONS
  % -----------------------------------------------------------------------
  % -----------------------------------------------------------------------
  
  % -----------------------------------------------------------------------
  %                       collect data
  % -----------------------------------------------------------------------
  d=[];
  tsl=[];
  drugID=[];
  nObs=nan(1,nDrugCondit);
  % loop over concs
  for dcIx=1:nDrugCondit
    % - start with detected parameters (in all variables: observations in
    % rows, parameters in columns)
    detD=cat(2,PSCR{g,drugIx(dcIx),detParIx});
    % - corresponding tsl
    detTsl=PSCR{g,drugIx(dcIx),allTslIx};
    % - get rid of entries with values outside the ones specified in fullPSCPar
    goodIx=true(size(detD,1),1);
    for ii=1:numel(detParIx)
      goodIx=goodIx & detD(:,ii)>=fullPSCPar{detParInd(ii),3}(1) & detD(:,ii)<fullPSCPar{detParInd(ii),3}(end);
    end
    detD=detD(goodIx,:);
    detTsl=detTsl(goodIx,:);
    
    % - now load fitted parameters
    fitD=cat(2,PSCR{g,drugIx(dcIx),fitParIx});
    % - corresponding tsl
    fitTsl=PSCR{g,drugIx(dcIx),tslIx};
    % check for nans etc
    goodIx=all(isfinite(fitD),2);
    fitD=fitD(goodIx,:);
    fitTsl=fitTsl(goodIx);
    
    % - now assign values of fitted PSC to proper entries in list of
    % detected PSCs:
    % -- identify PSCs by their time stamps
    [~,ia,ib]=intersect(detTsl,fitTsl);
    % -- preallocate new columns
    detD=cat(2,detD,nan(numel(detTsl),numel(fitParInd)));
    % -- assign
    detD(ia,numel(detParInd)+1:end)=fitD(ib,:);
    % - number of observations
    nObs(dcIx)=size(detD,1);
    % - drugID (=drug condition)
    curDrugID=repmat(drugIx(dcIx),[nObs(dcIx) 1]);
    % concatenate all
    d=cat(1,d,detD);
    drugID=cat(1,drugID,curDrugID);
    tsl=cat(1,tsl,detTsl);
  end
  
  % -----------------------------------------------------------------------
  %                compute clusters
  % -----------------------------------------------------------------------
  % general strategy:
  % - define clusters based on data from drug conditions listed in drugClusterIx
  % - assign data from other drug conditions to previously found clusters
  % - count events in each cluster separately for drug conditions
  % - assign bootstrapped data to previously found clusters
  % - count events
  
  % pick set of drug conditions to be compared statistically
  [~,tmpIx]=intersect(drugIx,drugStatsIx,'stable');
  % an array storing cluster IDs of each IPSC
  clustID=nan(sum(nObs),1);
  
  % - data to be used for clustering (all drug combinations)
  curD=d(:,plotPSCParInd);
  % - log transform?
  if doLogTransform
    curD=log10(curD);
  end
  % index to events in drug condition chosen for cluster generation
  curIx=ismember(drugID,drugClusterIx);

  switch clustMeth
    case 'kmeans'
      % evaluate number of clusters needed to partition the data set
      eva = evalclusters(curD(curIx,:),'kmeans','CalinskiHarabasz','KList',1:nMaxClust);
      disp(['optimal number of clusters: ' int2str(eva.OptimalK)]);
      nClust=nFixClust;
      [cid,clustCtr]=kmeans(curD(curIx,:),nClust,'distance',distType,'replicates',nReplic,'options',cluOpts);
      % store cluster IDs
      clustID(curIx)=cid;
      % now assign events in other drug conditions to clusters found (unless
      % clusters were determined based on all drug conditions to be dealt
      % with, in which case the job is done already)
      if ~isequal(drugClusterIx,drugIx)
        cid=kmeans(curD(~curIx,:),nClust,'distance',distType,'options',cluOpts,...
          'MaxIter',1,'Start',clustCtr);
        % store
        clustID(~curIx)=cid;
      end
      
    case 'hGauss'
      % eva = evalclusters(curD(curIx,:),'gmdistribution','CalinskiHarabasz','KList',1:nClust)
      gm = fitgmdist(curD(curIx,:),nClust,'Options',cluOpts);
      cid=cluster(gm,curD(curIx,:));
      clustCtr=gm.mu;

      % store cluster IDs
      clustID(curIx)=cid;
      
      % now assign events in other drug conditions to clusters found (unless
      % clusters were determined based on all drug conditions to be dealt
      % with, in which case the job is done already)
      if ~isequal(drugClusterIx,drugIx)
        cid=cluster(gm,curD(~curIx,:));
        % store
        clustID(~curIx)=cid;
      end
    otherwise
      error('illegal clustering method');
  end
  % §§ reorder clusters according to centroid values?
  
  
  % event counts per cluster:  cluster | drug condition 
  % ** note that numerical values in curDrugID are important
  clustCount=accumarray([clustID drugID],1);
  
  if doBootStats
    bootClustIx=[];
    % loop over drug conditions to be tested
    for dcIx=1:numel(drugStatsIx)
      curIx=find(drugID==drugStatsIx(dcIx));
      curD=d(curIx,plotPSCParInd);
      % - log transform?
      if doLogTransform
        curD=log(curD);
      end
      
      curNObs=numel(curIx);
      bootCurD=datasample(curD,curNObs*nBoot);
      % kmeans with iter=1 to assign bootstrapped data to clusters
      tmp=kmeans(bootCurD,nClust,'distance',distType,'options',cluOpts,...
        'MaxIter',1,'Start',clustCtr);
      % glue drugID array as second column
      tmp=cat(2,tmp,repmat(drugStatsIx(dcIx),size(tmp)));
      
      % glue array indicating bootstrapping instance
      tmp=cat(2,tmp,reshape(ones(curNObs,1)*(1:nBoot),[curNObs*nBoot 1]));
      
      % finally, glue things on top of each other
      bootClustIx=cat(1,bootClustIx,tmp);
    end
    bootClustCount=accumarray(bootClustIx,1);
    
    % keep columns to be analyzed and reshape for stats: boot instances | cluster | drug cond
    bootClustCount=permute(bootClustCount(:,drugStatsIx,:),[3 1 2]);
    % array indicating clusters
    tmp=repmat(1:nClust,[nBoot 1 numel(drugStatsIx)]);
    % array indicating drug conditions
    tmp2=repmat(permute(drugStatsIx,[1 3 2]),[nBoot nClust 1]);
    
    stats=mes2way(bootClustCount(:),[tmp(:) tmp2(:)],{'eta2'},...
      'fName',{'cluster','ACh condit'},'isDep',[1 0],'doDataPlot',false);
    
  end
  
  
  % -----------------------------------------------------------------------
  %  inhibitory impact, interpolated from events in that cluster which could be fitted
  % -----------------------------------------------------------------------
  helpInd=strcmp(fullPSCPar(:,1),'allAmp');
  ampSumArr=nan(nClust,nDrugCondit);
  for dcIx=1:nDrugCondit
    for cli=1:nClust
      % events of current drug condition and cluster
      curIx=clustID==cli & ismember(drugID,drugIx(dcIx));
      % amplitudes
      curD=d(curIx,helpInd);
      ampSumArr(cli,dcIx)=nansum(curD);
    end      
    %     % normalize: divide by total amplitude count
    %     ampSumArr(:,dcIx)=ampSumArr(:,dcIx)/sum(ampSumArr(:,dcIx));
  end
  
  % -----------------------------------------------------------------------
  %              assemble data & stats
  % -----------------------------------------------------------------------
  % assemble all interesting characteristics of clusters from ALL drug
  % conditions in one 4D variable:
  % cluster ID | drug condition | characteristics as detailed below |
  % experiment
  % - cluster centroid x coordinate (=mean of analysis par1 of that cluster)
  % - cluster centroid y coordinate (=mean of analysis par2 of that cluster)
  % - event counts in all drug conditions (one column per condition)
  % - 'inhibitory impact' in all drug conditions, quantified as the sum of
  % PSC amplitudes in cluster
  
  clustPar(:,:,:,g)=cat(3,...
    repmat(permute(clustCtr,[1 3 2]),[1 max(drugIx) 1]),... % note that this defines cluster centers as being the same across all drug conditions
    clustCount,...
    ampSumArr);
  
  % slice header
  clustParHeader=cat(2,{plotPSCPar{:,1}},{'PSC count','cumul. PSC impact'});
  
  % assemble all interesting omnibus analysis results:
  tabStats=mestab(clustCount(:,drugStatsIx));
  omniRes(g,1:3)=[tabStats.cramerV tabStats.cramerVCi];
  % column header
  omniResHeader={'V','ci(V) lo','ci(V) hi'};
  
  
  % -----------------------------------------------------------------------  
  %                       2D histograms
  % -----------------------------------------------------------------------  
  % each element of histPony contains the 2D histograms of the data
  % bins par1 | bins par2 | drug cond
  histPony=histTemplate;
  
  % - loop over concs
  for dcIx=1:nDrugCondit
    curIx=drugID==drugIx(dcIx);
    % ** retain only two columns as these are under scrutiny here
    curD=d(curIx,plotPSCParInd);
    curNObs=nObs(dcIx);
    
    if doHistScale
      hScaleFac=1/curNObs;
    else
      hScaleFac=1;
    end
    
    tmp=numel(find(curD(:,1)>plotPSCPar{1,3}(end)))/curNObs;
    if tmp>.01
      warning([num2str(tmp*100) ' % of events not covered by bin edges of par ' plotPSCPar{1,1}]);
    end
    tmp=numel(find(curD(:,2)>plotPSCPar{2,3}(end)))/curNObs;
    if tmp>.01
      warning([num2str(tmp*100) ' % of events not covered by bin edges of par ' plotPSCPar{2,1}]);
    end
    
    % § 2D histogram
    tmpH=hist3(curD,'edges',plotPSCPar(:,3))*hScaleFac;
    % tmpH=weighHist2D(curD,weight,plotPSCPar{c(pci,1),3},plotPSCPar{c(pci,2),3})*hScaleFac;
    
    % smooth by convolution with 2D Gaussian:
    tmpH=conv2(tmpH,gauss2d,'same');
    histPony(:,:,dcIx)=tmpH;
    
  end
  % compute 2D correlation as indicator of similarity
  [~,tmpIx]=intersect(drugIx,drugStatsIx,'stable');
  cc=corr2(histPony(:,:,tmpIx(1)),histPony(:,:,tmpIx(2)));
  % add to omniRes
  omniRes(g,4)=cc;
  % column header
  omniResHeader{4}='XC';
  
  
  % -----------------------------------------------------------------------
  % -----------------------------------------------------------------------
  %                       PLOTTING
  % -----------------------------------------------------------------------
  % -----------------------------------------------------------------------
  % ensure that plots are not too elongated
  nRow=3;
  nCol=max(3,nDrugCondit);
  
  mima=cellfun(@(x) x([1 end])',plotPSCPar(:,3),'UniformOutput',false);
  mima=cat(2,mima{:});
  
  figName=[LISTEXP{g} '_' plotPSCPar{1,1} '_' plotPSCPar{2,1}];
  fh=mkfig(g,'b');
  set(fh,'name',figName);
  colormap(coma('CubicL'));
  clf
  
  % --- SCATTER PLOTS ---
  subplot(nRow,nCol,1)
  hold on
  % raw points (plot for all drug conditions so that averages won't be overplotted)
  for dcIx=1:nDrugCondit
    curIx=drugID==drugIx(dcIx);
    ph=plot(d(curIx,1),d(curIx,2),'.');
    set(ph,'color',pCol(dcIx,:));
  end
  % - loop over drug conditions: averages & errors
  for dcIx=1:nDrugCondit
    curIx=drugID==drugIx(dcIx);
    if 0
      prc=prctile(d(curIx,:),[5 50 95]);
      av=prc(2,:);
      er=prc([1 3],:);
    else
      av=nanmean(d(curIx,:));
      st=nanstd(d(curIx,:));;
      er=[av-st; av+st];
    end
    
    lh=errorcross2(...
      [av(1),av(2)],...
      [av(1)-er(1,1)  av(2)-er(1,2)],...
      [er(2,1)-av(1)  er(2,2)-av(2)],...
      'color',pCol(dcIx,:)*.7,'linewidth',2.0);
    
    ph=plot(av(1),av(2),'o');
    set(ph,'color',pCol(dcIx,:)*.7,'markerfacecolor',pCol(dcIx,:),'markersize',8);
    
    if dcIx==nDrugCondit
      set(gca,'xscale','log','yscale','log');
      axis([mima(:,1)' mima(:,2)'])
      set(gca,'xtick',plotPSCPar{1,4},'ytick',plotPSCPar{2,4})
      grid on
      xlabel(axLabel{1});
      ylabel(axLabel{2});
    end
  end
  
  % show clusters in scatter plot
  subplot(nRow,nCol,2)
  hold on
  for cli=1:nClust
    % pick events only from drug condition chosen for cluster generation
    curIx=clustID==cli & ismember(drugID,drugClusterIx);
    ph=plot(d(curIx,1),d(curIx,2),'.');
    % set(ph,'color',clustCol(cli,:))
    th=text(10.^(clustCtr(cli,1)),10.^(clustCtr(cli,2)),int2str(cli),...
      'HorizontalAlignment', 'center','fontsize',15,'fontweight','b');
  end
  set(gca,'xscale','log','yscale','log');
  axis([mima(:,1)' mima(:,2)'])
  set(gca,'xtick',plotPSCPar{1,4},'ytick',plotPSCPar{2,4})
  grid on
  xlabel(axLabel{1});
  ylabel(axLabel{2});
  box off
  
 % bubble plot - amplitudes
  subplot(nRow,nCol,3)
  hold on
  for dcIx=1:nDrugCondit
    curIx=drugID==drugIx(dcIx);
    ph=scatter(d(curIx,1),d(curIx,2),d(curIx,strcmp(fullPSCPar(:,1),'allAmp'))*ampScalFac);
    set(ph,'MarkerEdgeColor',pCol(dcIx,:));
  end
  set(gca,'xscale','log','yscale','log');
  axis([mima(:,1)' mima(:,2)'])
  set(gca,'xtick',plotPSCPar{1,4},'ytick',plotPSCPar{2,4})
  grid on
  xlabel(axLabel{1});
  ylabel(axLabel{2});
  box off
  
  % --- SURFACE PLOTS ---
  if doSameColorScale
    clim=prctile(histPony(:),[.5 99.9]);
  end
  for dcIx=1:nDrugCondit
    sph=subplot(nRow,nCol,nCol+dcIx);
    if any(any(isfinite(histPony(:,:,dcIx))))
      % unfortunately, grid does not show up in surface plot viewed from
      % above
      % h=surf(plotPSCPar{1,3},plotPSCPar{2,3},histPony(:,:,dcIx)');
      % view(0,90);
      [~,h]=contourf(plotPSCPar{1,3},plotPSCPar{2,3},histPony(:,:,dcIx)',40);
      set(h,'linestyle','none');
      if doSameColorScale
        set(gca,'cLim',clim);
      end
      set(gca,'xscale','log','yscale','log');
      set(gca,'xtick',plotPSCPar{1,4},'ytick',plotPSCPar{2,4})
      grid on
      grid minor
      sph.GridAlpha=.8;
      sph.GridLineStyle='--';
      % same axis limits as corresponding scatter plot
      axis([mima(:,1)' mima(:,2)'])
      % info on IPSC freq
      smarttext(['f=' num2str(PSCR{g,drugIx(dcIx),freqParIx(1)},2) '/'...
        num2str(PSCR{g,drugIx(dcIx),freqParIx(2)},2) ' Hz'],.0,.8,'color','w','fontsize',10,'fontweight','bold');
      xlabel(axLabel{1});
    end
    title(drugTag{dcIx},'fontsize',12);
  end
  smarttext(['cc=' num2str(omniRes(g,strcmp(omniResHeader,'XC')),2)],.0,.5,'color','w','fontsize',10,'fontweight','bold');
  
  % difference of UNSCALED histograms, divided by their sum
  sph=subplot(nRow,nCol,2*nCol+1);
  colormap(sph,coma('parula'));
  [~,tmpIx]=intersect(drugIx,drugStatsIx,'stable');
  A=histPony(:,:,tmpIx(1))*nObs(tmpIx(1));
  B=histPony(:,:,tmpIx(2))*nObs(tmpIx(2));
  tmpH=(A-B)./(A+B);
  % get rid of nans and infs
  tmpH(~isfinite(tmpH))=0;
  if any(any(isfinite(tmpH)))
    % symmetric color limits
    % cl=max(abs(tmpH(:)))*[-1 1];
    cl=[-1 1];
    [~,h]=contourf(plotPSCPar{1,3},plotPSCPar{2,3},tmpH',40);
    set(h,'linestyle','none');
    set(gca,'clim',cl);
    set(gca,'xscale','log','yscale','log');
    set(gca,'xtick',plotPSCPar{1,4},'ytick',plotPSCPar{2,4})
    grid on
    grid minor
    sph.GridAlpha=.8;
    sph.GridLineStyle='--';
    % same axis limits as corresponding scatter plot
    axis([mima(:,1)' mima(:,2)'])
    
    title('diff')
  end

  % plot number of events in each cluster and display Cramer's V
  sph=subplot(nRow,nCol,2*nCol+2);
  hold on
  for dcIx=1:numel(drugStatsIx)
    ph=plot(clustCount(:,drugStatsIx(dcIx)),'k-o');
    set(ph,'markerfacecolor',pCol(drugStatsIx(dcIx),:),'markersize',8);
  end
  nicexy0ax(8);
  % logarithmic y axis so we can judge proportions
  set(gca,'xtick',1:nClust,'yscale','log')
  smarttext(['Cramer''s V=' num2str(tabStats.cramerV,2)],.0,.8,'color','k','fontweight','bold');
  grid on
  xlabel('cluster number')
  ylabel('event counts');

 % bubble plot of clusters - inhibitory amplitudes 
  subplot(nRow,nCol,9)
  hold on
  for dcIx=1:nDrugCondit
    ph=scatter(10.^clustCtr(:,1),10.^clustCtr(:,2),ampSumArr(:,dcIx)*ampScalFac/10,'filled');
    set(ph,'MarkerEdgeColor',pCol(dcIx,:)*.5,'MarkerFaceColor',pCol(dcIx,:));
  end
  set(gca,'xscale','log','yscale','log');
  axis([mima(:,1)' mima(:,2)'])
  set(gca,'xtick',plotPSCPar{1,4},'ytick',plotPSCPar{2,4})
  grid on
  xlabel(axLabel{1});
  ylabel(axLabel{2});
  box off

  
  
%   subplot(nRow,nCol,9)
%   hold on
%   helpInd=strcmp(fullPSCPar(:,1),'tDecay');
%   for cli=1:nClust
%     % pick events only from drug condition chosen for cluster generation
%     curIx=clustID==cli & ismember(drugID,drugClusterIx);
%     n=histcounts(d(curIx,helpInd),fullPSCPar{helpInd,3});
%     ph=stairs(fullPSCPar{helpInd,3}(1:end-1),n);
%     % set(ph,'color',clustCol(cli,:))
%   end
%   set(gca,'xscale','log','yscale','log');
%   grid on
%   xlabel('tau');
%   ylabel('N');
%   
  
  drawnow
  pause(.05)

  % ----- printing -----
  if ~isempty(printas),
    print([plotPath '\' figName],printas);
  end
  
end

% % dump all results of interest in base workspace
% namePart=dFn(1:findstr(dFn,'.')-1);
% assignin('base',['clustPar_' namePart],clustPar);
% assignin('base','clustParHeader',clustParHeader);
% assignin('base',['omniRes_' namePart],omniRes);
% assignin('base','omniResHeader',omniResHeader);

namePart=dFn(1:findstr(dFn,'.')-1);
save([dataPath dSubDir namePart '_clust'],'clustPar*','omniRes*');
% =========================================================================
% --- LOCAL FUNCTIONS --- LOCAL FUNCTIONS --- LOCAL FUNCTIONS --- LOCAL FUN
% =========================================================================

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
    