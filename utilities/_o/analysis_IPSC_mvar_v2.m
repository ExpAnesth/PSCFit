function analysis_IPSC_mvar(ds)
% ** function analysis_IPSC_mvar(ds)
% multivariate analysis of IPSCs 


% to do:
% MAJOR
% - maybe implement plotting raw IPSCs in selected bin(s)?
% SMALL FRY



% -------------------- PREPARATIONS ---------------------------------------

doBootStats=ds.nBoot>0;
anPSCPar=ds.fullPSCPar(ds.anPSCParInd,:);

% -------------------------------------------------------------------------
% --------------------- run, rabbit, run! ---------------------------------
% -------------------------------------------------------------------------

% precompute some frequently needed variables
nDrugCondit=numel(ds.drugIx);
nFullPSCPar=size(ds.fullPSCPar,1);
nAnPSCPar=size(anPSCPar,1);
axLabel=anPSCPar(1:nAnPSCPar,2);
% number of nonredundant combinations of PCs/parameters to analyze and
% plot
[c,nc]=combin(nAnPSCPar,'autoC',0);
% 2d Gaussian for smoothing
gauss2d=Gauss2D;

% LOAD DATA (** be aware that the files also contain struct ds **)
load([ds.dataPath ds.dSubDir ds.dFn],'PSCR','PSCRMN','depPar','expDate','LISTEXP');
nExp=size(PSCR,1);
% indexes into data:
% - analyzed parameters (=slices of PSCR)
[~,~,parIx]=intersect(anPSCPar(:,1)',depPar,'stable');
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
detParInd=find(strcmp(ds.fullPSCPar(:,5),'detected'));
fitParInd=find(strcmp(ds.fullPSCPar(:,5),'fitted'));
[~,~,detParIx]=intersect(ds.fullPSCPar(detParInd,1)',depPar,'stable');
[~,~,fitParIx]=intersect(ds.fullPSCPar(fitParInd,1)',depPar,'stable');

% preallocate table of 'omnibus' results
% - construct column names: 'r_' for correlation coeff, then parameter name
% (e.g. tau), then either nothing or 'ciLo' or 'ciUp' for lower and upper
% confidence intervals
tmpVn=cellfun(@(c) cat(2,'r_',c),ds.fullPSCPar(ds.anPSCParInd2,1)','uniformoutput',false);
tmpVn(end+1,:)=cellfun(@(c) cat(2,c,'_ciLo'),tmpVn(1,:),'uniformoutput',false);
tmpVn(end+1,:)=cellfun(@(c) cat(2,c,'_ciUp'),tmpVn(1,:),'uniformoutput',false);

tmpVn1=cellfun(@(c) cat(2,'loPrc_',c),ds.fullPSCPar(ds.anPSCParInd2,1)','uniformoutput',false);
tmpVn2=cellfun(@(c) cat(2,'mdPrc_',c),ds.fullPSCPar(ds.anPSCParInd2,1)','uniformoutput',false);
tmpVn3=cellfun(@(c) cat(2,'hiPrc_',c),ds.fullPSCPar(ds.anPSCParInd2,1)','uniformoutput',false);

varNm={'distXC','distXC_z','distB','distB_z','loPrc_prop','mdPrc_prop','hiPrc_prop',...
  tmpVn1{:},tmpVn2{:},tmpVn3{:},tmpVn{:}};

nanCell=cell(1,numel(varNm));
[nanCell{:}]=deal(nan(nExp,1));
omniRes=table(nanCell{:},'variablenames',varNm);

% loop over experiments
for g=1:nExp
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
    detD=cat(2,PSCR{g,ds.drugIx(dcIx),detParIx});
    % - corresponding tsl
    detTsl=PSCR{g,ds.drugIx(dcIx),allTslIx};
    % - get rid of entries with values outside the ones specified in ds.fullPSCPar
    goodIx=true(size(detD,1),1);
    for ii=1:numel(detParIx)
      goodIx=goodIx & detD(:,ii)>=ds.fullPSCPar{detParInd(ii),3}(1) & detD(:,ii)<ds.fullPSCPar{detParInd(ii),3}(end);
    end
    detD=detD(goodIx,:);
    detTsl=detTsl(goodIx,:);
    
    % - now load fitted parameters
    fitD=cat(2,PSCR{g,ds.drugIx(dcIx),fitParIx});
    % - corresponding tsl
    fitTsl=PSCR{g,ds.drugIx(dcIx),tslIx};
    % check for nans etc
    goodIx=all(isfinite(fitD),2);
    fitD=fitD(goodIx,:);
    fitTsl=fitTsl(goodIx);
    
    % *** now assign values of fitted PSC to proper entries in list of
    % detected PSCs ***
    % -- identify PSCs by their time stamps
    [~,ia,ib]=intersect(detTsl,fitTsl);
    % -- preallocate new columns
    detD=cat(2,detD,nan(numel(detTsl),numel(fitParInd)));
    % -- assign
    detD(ia,numel(detParInd)+1:end)=fitD(ib,:);
    
    % loop over parameters, kicking out values outside of bounds defined by
    % bins in anPSCPar (copied from ds.fullPSCPar):
    % - total number of observations (determined by number of 'detected'
    % events)
    tmpNObs=size(detD,1);
    % - number of valid observations with respect to two parameters to be
    % analyzed (which may differ from the above if a 'fitted' parameter is
    % included)
    tmpNGoodObs=sum(all(isfinite(detD(:,ds.anPSCParInd)),2));
    goodIx=true(tmpNObs,1);
    for ii=1:nAnPSCPar
      goodIx=goodIx & ...
        detD(:,ds.anPSCParInd(ii))>=anPSCPar{ii,3}(1) &...
        detD(:,ds.anPSCParInd(ii))<=anPSCPar{ii,3}(end);
    end
    disp([num2str(sum(goodIx)/sum(tmpNGoodObs)*100,3) ' % of events within bounds of bins'])
    detD=detD(goodIx,:);
    % - number of observations
    nObs(dcIx)=sum(goodIx);
    % - drugID (=drug condition)
    curDrugID=repmat(ds.drugIx(dcIx),[nObs(dcIx) 1]);
    % concatenate all
    d=cat(1,d,detD);
    drugID=cat(1,drugID,curDrugID);
    tsl=cat(1,tsl,detTsl);
  end
  nObsAll=sum(nObs);
  % set up array containing indices to bins to which individual IPSCs
  % belong
  binID=nan(nObsAll,2);
  
  % if bootstrapping or random permutation is requested...
  if doBootStats
    % preallocate
    d(:,:,2:ds.nBoot+1)=nan;
    if 0
      % lazy way: draw out of total population
      for bIx=1:ds.nBoot
        d(:,:,bIx+1)=d(randperm(nObsAll),:,1);
      end
    else
      % more correct way:
      % make temporary subset of d, comprised of populations to be compared
      curIx=ismember(drugID,ds.drugStatsIx);
      tmpD=d(curIx,:,1);
      tmpNObs=sum(nObs(ds.drugStatsIx));
      % fill value for population(s) not to be compared (so that accumarray
      % further below runs)
      filla=repmat(nanmean(tmpD),[sum(~curIx) 1]);
      for bIx=1:ds.nBoot
        d(curIx,:,bIx+1)=tmpD(randperm(tmpNObs),:);
        d(~curIx,:,bIx+1)=filla;
      end
    end
    clear tmp* filla
  end
  
  % -----------------------------------------------------------------------  
  %         computation of 2D histogram & derived parameters
  % -----------------------------------------------------------------------  
  
  % 2D histogram: par1 bin edges | par2 bin edges | original & bootstrapped | drug conditions
  histArr=nan(numel(anPSCPar{1,3})-1,numel(anPSCPar{2,3})-1,ds.nBoot+1,nDrugCondit);

  % - loop over concs
  for dcIx=1:nDrugCondit
    % pick events of current drug condition
    curIx=drugID==ds.drugIx(dcIx);
    % ** retain only two columns as these are under scrutiny here
    curD=d(curIx,ds.anPSCParInd,:);
    curNObs=nObs(dcIx);
    % check bin coverage
    for ii=1:2
      tmp=numel(find( curD(:,ii,1)>anPSCPar{ii,3}(end) | curD(:,ii,1)<anPSCPar{ii,3}(1) ))/curNObs;
      if tmp>.01
        warning([num2str(tmp*100) ' % of events not covered by bin edges of par ' anPSCPar{ii,1}]);
      end
    end
       
%     tmpH=weighHist2D(curD,weight,anPSCPar{c(pci,1),3},anPSCPar{c(pci,2),3})*hScaleFac;
    
    % use functions discretize and accumarray for fast construction of 2D
    % histograms of original and shuffled data
    for ii=1:2
      curD(:,ii,:)=discretize(curD(:,ii,:),anPSCPar{ii,3});
      % store bin indexes
      binID(curIx,ii)=curD(:,ii,1);
    end
    % helper array coding for third dim in histArr, namely, the bootstrap
    % instance
    tmpArr=ones(curNObs,1)*(1:ds.nBoot+1);
    % in order to get the 'slices' of curD stacked in the first dim via
    % reshape we first have to exchange 2nd and 3rd dim
    curD=permute(curD,[1 3 2]);
    histArr(:,:,:,dcIx)=accumarray([reshape(curD,[curNObs*(ds.nBoot+1) 2]) tmpArr(:)],...
      1,size(histArr(:,:,:,dcIx)));
    % smooth by convolution with 2D Gaussian (alas, no alternative to a
    % loop here)
    for ii=1:ds.nBoot+1
      histArr(:,:,ii,dcIx)=conv2(histArr(:,:,ii,dcIx),gauss2d,'same');
    end    
  end
  
  % compute various indexes of similarity/distance:
  % - loop over [original + bootstrap instances] for stats
  cc=nan(ds.nBoot+1,1);
  bd=nan(ds.nBoot+1,1);
  for ii=1:ds.nBoot+1
    % extract and normalize (** we have to normalize by the individual sum
    % due to Gaussian smoothing)
    tmpD1=histArr(:,:,ii,ds.drugStatsIx(1));
    tmpD1=tmpD1/sum(tmpD1(:));
    tmpD2=histArr(:,:,ii,ds.drugStatsIx(2));
    tmpD2=tmpD2/sum(tmpD2(:));
    % - 2D correlation
    cc(ii)=corr2(tmpD1,tmpD2);
    % - Bhattacharayya distance
    bd(ii)=sqrt(1-sum(sqrt(tmpD1(:).*tmpD2(:))));
  end
  % add value of comparison of originals to omniRes
  [~,ia]=intersect(omniRes.Properties.VariableNames,{'distXC'});
  omniRes(g,ia)={cc(1)};
  % compute z score
  [~,mu,sig]=zscore(cc(2:end));
  % add z score
  [~,ia]=intersect(omniRes.Properties.VariableNames,{'distXC_z'});
  omniRes(g,ia)={(cc(1)-mu)/sig};
  
  % add value of comparison of originals to omniRes  
  [~,ia]=intersect(omniRes.Properties.VariableNames,{'distB'});
  omniRes(g,ia)={bd(1)};
  % compute z score
  [~,mu,sig]=zscore(bd(2:end));
  % add z score
  [~,ia]=intersect(omniRes.Properties.VariableNames,{'distB_z'});
  omniRes(g,ia)={(bd(1)-mu)/sig};
  
  % compute bin-wise normalized difference of histograms: (A-B)/(A+B)
  bwNormHistDiff=(histArr(:,:,:,ds.drugStatsIx(1))-histArr(:,:,:,ds.drugStatsIx(2)))./sum(histArr(:,:,:,ds.drugStatsIx),4);
  % exclude bins which contain less than threshNbinwise events by setting
  % their value to neutral (zero); as the histograms are smoothed, set it to
  % a value of [minimal required number]*[peak of Gaussian]
  threshNbinwise=max(gauss2d(:))*1;
  badIx=sum(histArr(:,:,1,ds.drugStatsIx),4)<=threshNbinwise;
  bwNormHistDiff(repmat(badIx,[1 1 ds.nBoot+1]))=nan;

  % convert binID (identity of bin in which each IPSC resides) from subindex to linear index
  binID=sub2ind([numel(anPSCPar{1,3})-1,numel(anPSCPar{2,3})-1],binID(:,1),binID(:,2));
  
  % compute at which percentile the original value is of the shuffled
  % distribution
  bwHistDiffIndex=(sum(bwNormHistDiff(:,:,2:end)<repmat(bwNormHistDiff(:,:,1),[1 1 ds.nBoot]),3) + ...
    0.5*sum(bwNormHistDiff(:,:,2:end)==repmat(bwNormHistDiff(:,:,1),[1 1 ds.nBoot]),3))/ds.nBoot*100;
  
  % § set to value which is neutral in current version of plots
  bwHistDiffIndex(badIx)=50;
  
  % assign to each PSC a 'difference score' which is simply the value of
  % bwNormHistDiff of the bin in which the PSC belongs (the higher the
  % value, the more are PSCs in the first drug condition (as listed in
  % ds.drugStatsIx) represented in that bin) (note linear indexing of a 3D
  % variable)
  diffScore=bwNormHistDiff(binID);
  % add percentile of null distribution to which it belongs
  diffScore(:,2)=bwHistDiffIndex(binID);
  
  
  % *** compute properties of IPSCs in different groups
  % identify IPSCs in 'depressed' and 'induced' bins in first condition
  % listed in ds.drugStatsIx
  curIx=drugID==ds.drugStatsIx(1);

  for ii=1:3
    switch ii
      case 1
        % lower percentiles
        str='loPrc_';
        tmpIx=diffScore(:,2)<=ds.ndPrctile{1,1};
      case 2
        % upper percentiles
        tmpIx=diffScore(:,2)>=ds.ndPrctile{2,1};
        str='hiPrc_';
      case 3
        % all in between
        tmpIx=diffScore(:,2)>ds.ndPrctile{1,1} & diffScore(:,2)<ds.ndPrctile{2,1};
        str='mdPrc_';
    end
    
    % proportion of IPSCs
    [~,ia]=intersect(omniRes.Properties.VariableNames,{[str 'prop']});
    omniRes(g,ia)={sum(curIx & tmpIx)/nObs(ds.drugStatsIx(1))};
    
    for pIx=1:numel(ds.anPSCParInd2)
      [~,ia]=intersect(omniRes.Properties.VariableNames,{[str ds.fullPSCPar{ds.anPSCParInd2(pIx),1}]});
      tmpD=nanmedian(d(curIx & tmpIx,ds.anPSCParInd2(pIx)));
      omniRes(g,ia)={tmpD};
    end
  end
  
  % compute correlation between 'difference score' and log of parameters
  % specified in ds.anPSCParInd2
  for pIx=1:numel(ds.anPSCParInd2)
    [r,~,rl,ru]=corrcoef(diffScore(curIx,2),log(d(curIx,ds.anPSCParInd2(pIx))),'rows','complete');
    [~,ia]=intersect(omniRes.Properties.VariableNames,{['r_'   ds.fullPSCPar{ds.anPSCParInd2(pIx),1}]});    
    omniRes(g,ia)={r(1,2)};
    [~,ia]=intersect(omniRes.Properties.VariableNames,{['r_'   ds.fullPSCPar{ds.anPSCParInd2(pIx),1} '_ciLo']});    
    omniRes(g,ia)={rl(1,2)};
    [~,ia]=intersect(omniRes.Properties.VariableNames,{['r_'   ds.fullPSCPar{ds.anPSCParInd2(pIx),1} '_ciUp']});    
    omniRes(g,ia)={ru(1,2)};
  end
  
  % -----------------------------------------------------------------------
  % -----------------------------------------------------------------------
  %                       PLOTTING
  % -----------------------------------------------------------------------
  % -----------------------------------------------------------------------
  
  if ds.doPlot
    
    % ensure that plots are not too elongated
    nRow=4;
    nCol=4;
    % axes for 2D histogram plots
    xax=conv(anPSCPar{1,3},[1 1]*.5,'valid');
    yax=conv(anPSCPar{2,3},[1 1]*.5,'valid');
    
    % we need a little helper index
    pInd=ds.anPSCParInd;
    
    % axis limits
    if ds.doIndividualAxLim
      mima=[min(d(:,pInd)); max(d(:,pInd))];
    else
      mima=cellfun(@(x) x([1 end])',anPSCPar(:,3),'UniformOutput',false);
      mima=cat(2,mima{:});
    end
    
    figName=[LISTEXP{g} '_' anPSCPar{1,1} '_' anPSCPar{2,1}];
    fh=mkfig(g,'b');
    orient landscape
    set(fh,'name',figName);
    colormap(coma('CubicL'));
    clf
    
    
    % ------------------------ SCATTER PLOTS --------------------------------
    % - scatter plot with error crosses
    subplot(nRow,nCol,1)
    hold on
    phArr=gobjects(1,nDrugCondit);
    % raw points (plot for all drug conditions so that averages won't be overplotted)
    for dcIx=1:nDrugCondit
      curIx=drugID==ds.drugIx(dcIx);
      ph=plot(d(curIx,pInd(1)),d(curIx,pInd(2)),'.');
      set(ph,'color',ds.pCol(dcIx,:));
    end
    % - loop over drug conditions: averages & errors
    for dcIx=1:nDrugCondit
      curIx=drugID==ds.drugIx(dcIx);
      if 1
        prc=prctile(d(curIx,pInd),[25 50 75]);
        av=prc(2,:);
        er=prc([1 3],:);
      else
        av=nanmean(d(curIx,pInd));
        st=nanstd(d(curIx,pInd));
        er=[av-st; av+st];
      end
      
      lh=errorcross2(...
        [av(1),av(2)],...
        [av(1)-er(1,1)  av(2)-er(1,2)],...
        [er(2,1)-av(1)  er(2,2)-av(2)],...
        'color',ds.pCol(dcIx,:)*.7,'linewidth',2.0);
      
      ph=plot(av(1),av(2),'o');
      set(ph,'color',ds.pCol(dcIx,:)*.7,'markerfacecolor',ds.pCol(dcIx,:),'markersize',6);
      phArr(dcIx)=ph(1);
      
      if dcIx==nDrugCondit
        set(gca,'xscale','log','yscale','log');
        axis([mima(:,1)' mima(:,2)'])
        set(gca,'xtick',anPSCPar{1,4},'ytick',anPSCPar{2,4})
        grid on
        xlabel(axLabel{1});
        ylabel(axLabel{2});
        lh=legend(phArr,ds.drugTag);
      end
    end
    colorbar
    
    % ---------------- SURFACE PLOTS ----------------------------------------
    if ds.doSameColorScale
      tmpHist=histArr(:,:,1,:);
      clim=[0 prctile(tmpHist(:),99.9)];
    end
    
    % histograms at all analyzed drug conditions
    for dcIx=1:nDrugCondit
      sph=subplot(nRow,nCol,nCol+dcIx);
      if any(any(isfinite(histArr(:,:,1,dcIx))))
        % unfortunately, grid does not show up in surface plot viewed from
        % above
        % h=surf(xax,yax,histArr(:,:,dcIx)');
        % view(0,90);
        [~,h]=contourf(xax,yax,histArr(:,:,1,dcIx)',20);
        set(h,'linestyle','none');
        if ds.doSameColorScale
          set(gca,'cLim',clim);
        end
        set(gca,'xscale','log','yscale','log');
        set(gca,'xtick',anPSCPar{1,4},'ytick',anPSCPar{2,4})
        grid on
        grid minor
        sph.GridAlpha=.8;
        sph.GridLineStyle='--';
        % same axis limits as corresponding scatter plot
        axis([mima(:,1)' mima(:,2)'])
        % info on IPSC freq
        smarttext(['f=' num2str(PSCR{g,ds.drugIx(dcIx),freqParIx(1)},2) '/'...
          num2str(PSCR{g,ds.drugIx(dcIx),freqParIx(2)},2) ' Hz'],.0,.8,'color','w','fontsize',8,'fontweight','bold');
        xlabel(axLabel{1});
        ylabel(axLabel{2});
      end
      title(ds.drugTag{dcIx});
      cbh=colorbar('location','eastoutside');
    end
    %   % display value of Bhatt. distance?
    %   [~,ia]=intersect(omniRes.Properties.VariableNames,{'distB'});
    %   smarttext(['Bhatt. dist=' num2str(omniRes{g,ia},2)],.0,.5,'color','w','fontsize',10,'fontweight','bold');
    
    % normalized difference of histograms
    sph=subplot(nRow,nCol,2*nCol+1);
    colormap(sph,coma('bluered'));
    tmpH=bwNormHistDiff(:,:,1);
    if any(any(isfinite(tmpH)))
      % symmetric color limits
      cl=[-1 1];
      [~,h]=contourf(xax,yax,tmpH',20);
      set(h,'linestyle','none');
      % overplot contour lines of 'significant' differences
      hold on
      for ii=1:size(ds.ndPrctile,1)
        [~,h2]=contour(xax,yax,bwHistDiffIndex',ds.ndPrctile{ii,1}*[1 1]);
        set(h2,'linewidth',1,'color',ds.ndPrctile{ii,2});
      end
      set(gca,'clim',cl);
      set(gca,'xscale','log','yscale','log');
      set(gca,'xtick',anPSCPar{1,4},'ytick',anPSCPar{2,4})
      grid on
      grid minor
      sph.GridAlpha=.8;
      sph.GridLineStyle='--';
      colorbar('location','eastoutside');
      % same axis limits as corresponding scatter plot
      axis([mima(:,1)' mima(:,2)'])
      xlabel(axLabel{1});
      ylabel(axLabel{2});
      title('norm. difference')
    end
    
    % example of histogram resulting from random choice of union
    for dcIx=1 %:numel(ds.drugStatsIx);
      sph=subplot(nRow,nCol,2*nCol+2);
      if any(any(isfinite(histArr(:,:,1,ds.drugStatsIx(dcIx)))))
        [~,h]=contourf(xax,yax,histArr(:,:,2,ds.drugStatsIx(dcIx))',20);
        set(h,'linestyle','none');
        if ds.doSameColorScale
          set(gca,'cLim',clim);
        end
        % mark exemplary bin
        hold on
        % ph=plot(mean(xax(exBinIx(1)+[0 1])),mean(yax(exBinIx(1)+[0 1])),'w*');
        % set(ph,'markersize',8);
        lh=line(xax(ds.exBinIx(1)+[0 1])'*[1 1],repmat(yax(ds.exBinIx(2)+[0 1]),2,1),'color','w','linewidth',1.5);
        lh=line(repmat(xax(ds.exBinIx(1)+[0 1]),2,1),yax(ds.exBinIx(2)+[0 1])'*[1 1],'color','w','linewidth',1.5);
        set(gca,'xscale','log','yscale','log');
        set(gca,'xtick',anPSCPar{1,4},'ytick',anPSCPar{2,4})
        grid on
        grid minor
        sph.GridAlpha=.8;
        sph.GridLineStyle='--';
        % same axis limits as corresponding scatter plot
        axis([mima(:,1)' mima(:,2)'])
        xlabel(axLabel{1});
        ylabel(axLabel{2});
      end
      title('random');
      cbh=colorbar('location','eastoutside');
    end
    
    % for one bin, show distribution of normalized count difference as well
    % as expected value given IPSC frequencies and assuming that
    % distributions are preserved
    sph=subplot(nRow,nCol,2*nCol+3);
    hold on
    h=histogram(bwNormHistDiff(ds.exBinIx(1),ds.exBinIx(2),2:end),'BinMethod','sqrt',...
      'displaystyle','stairs','edgecolor','k');
    %   h.FaceColor='k';
    %   h.FaceAlpha=1;
    nicexy0ax;
    yl=get(gca,'ylim');
    ph=plot(bwNormHistDiff(ds.exBinIx(1),ds.exBinIx(2),1),yl(2)*.9,'kv');
    set(ph,'markerfacecolor','k','markersize',7);
    % compute and plot expected value
    a=PSCR{g,ds.drugStatsIx(1),freqParIx(1)};
    b=PSCR{g,ds.drugStatsIx(2),freqParIx(1)};
    expVal=(a-b)/(a+b);
    lh=line(expVal*[1 1],yl,'linewidth',1,'color','b');
    nicexy0ax;
    xlabel('norm. count diff.');
    ylabel('N')
    
    % ------------------------ OTHER PLOTS ----------------------------------
    % cumulative density function of IPSCs vs percentiles of difference
    % scores
    sph=subplot(nRow,nCol,2*nCol+4);
    hold on
    %   for ii=1:2
    %     lh=line(ds.ndPrctile{ii,1}*[1 1],[0 1],'linestyle','--','color',ds.ndPrctile{ii,2});
    %   end
    paH=patch([0 0 ds.ndPrctile{1,1}*[1 1]],[0 1 1 0],ds.ndPrctile{1,2});
    paH.EdgeColor=ds.ndPrctile{1,2};
    paH=patch([ds.ndPrctile{2,1}*[1 1] [1 1]*100],[0 1 1 0],ds.ndPrctile{2,2});
    paH.EdgeColor=ds.ndPrctile{2,2};
    
    % index to drug cond of interest
    curIx=drugID==ds.drugStatsIx(1);
    h=cdfplot(diffScore(curIx,2));     % ecdf(diffScore(curIx,2))
    title('')
    set(h,'linewidth',.75,'color','k');
    xlabel('percentile of null distrib.')
    ylabel('cum. prob.')
    
    
    % scatter plots of PSC parameters vs. 'difference scores'
    % - we need a helper index to point to the PSC parameters to plot
    pInd2=ds.anPSCParInd2;
    nPInd2=numel(pInd2);
    phArr=gobjects(nPInd2,1);
    for pIx=1:nPInd2
      sph=subplot(nRow,nCol,3*nCol+pIx);
      hold on
      for dcIx=1:1 %numel(ds.drugStatsIx)
        % index to drug cond of interest
        curIx=drugID==ds.drugStatsIx(dcIx);
        % plot all
        ph=plot(diffScore(curIx,2),d(curIx,pInd2(pIx)),'ko');
        set(ph,'markerfacecolor',[.9 .9 .9],'linewidth',.05);
        % overplot PSCs in bins with percentile below/above threshold
        for ii=1:2
          if ii==1
            tmpIx=diffScore(:,2)<=ds.ndPrctile{ii,1};
          else
            tmpIx=diffScore(:,2)>=ds.ndPrctile{ii,1};
          end
          ph=plot(diffScore(curIx & tmpIx,2),d(curIx & tmpIx,pInd2(pIx)),'ko');
          set(ph,'markerfacecolor',ds.ndPrctile{ii,2},'linewidth',.05);
        end
        %         % compute correlation (all)
        %         [r,~,rl,ru]=corrcoef(diffScore(curIx,1),log(d(curIx,pInd2(pIx))),'rows','complete');
      end
      set(gca,'yscale','log')
      axis tight
      % make sure symbols don't protrude
      set(gca,'xlim',[-2 102]);
      %       title(['r=' num2str(r(1,2),3) ' [' num2str([rl(1,2),ru(1,2)],3) ']']);
      % xlabel('norm. count diff.');
      xlabel('percentile of null distrib.')      
      ylabel(ds.fullPSCPar{pInd2(pIx),2});
    end
    
    drawnow
    pause(.05)
    
    % ----- printing -----
    if ~isempty(ds.printas),
      print([ds.plotPath '\' figName],'-painters','-r400',ds.printas); % 
    end
  end
end

namePart=ds.dFn(1:findstr(ds.dFn,'.')-1);
save([ds.dataPath ds.dSubDir namePart '_mvar' int2str(ds.drugStatsIx')'],'omniRes*','ds');

assignin('base','omniRes',omniRes)

% =========================================================================
% --- LOCAL FUNCTIONS --- LOCAL FUNCTIONS --- LOCAL FUNCTIONS --- LOCAL FUN
% =========================================================================

function g=Gauss2D
% support up to 3 sigma
mu = [0 0];
Sigma = [1 0;0 1]*0.5;
% dimensions of muetze
x1 = -1:1:1; 
x2 = -1:1:1;
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

% -------------------------------------------------------------------------
%                         outdated plots 
% -------------------------------------------------------------------------

  %   % - bubble plot - amplitudes
  %   subplot(nRow,nCol,3)
  %   hold on
  %   for dcIx=1:nDrugCondit
  %     curIx=drugID==ds.drugIx(dcIx);
  %     ph=scatter(d(curIx,pInd(1)),d(curIx,pInd(2)),d(curIx,strcmp(ds.fullPSCPar(:,1),'allAmp'))*ds.ampScalFac);
  %     set(ph,'MarkerEdgeColor',ds.pCol(dcIx,:));
  %   end
  %   set(gca,'xscale','log','yscale','log');
  %   axis([mima(:,1)' mima(:,2)'])
  %   set(gca,'xtick',anPSCPar{1,4},'ytick',anPSCPar{2,4})
  %   grid on
  %   xlabel(axLabel{1});
  %   ylabel(axLabel{2});
  %   box off

  

%     % show distribution of distance values
%   sph=subplot(nRow,nCol,2*nCol+dcIx+1);
%   hold on
%   h=histogram(bd(2:end),'displaystyle','stairs','edgecolor','k');
%   %   h.FaceColor='k';
%   %   h.FaceAlpha=1;
%   ph=plot(bd(1),5,'k^');
%   set(ph,'markerfacecolor','w','markersize',8);
%   nicexy0ax;
%   xlabel('distance ctrl - ACh^+');
%   ylabel('N')
  
    %  % - bubble plot - decay times
    %   subplot(nRow,nCol,3)
    %   hold on
    %   % colormap coding for decay time
    %   cm=coma('turquoisamber','n',128);
    %   helpInd=strcmp(ds.fullPSCPar(:,1),'tDecay');
    %   % map decay times to colormap, clipping extremes
    %   maxDecay=20;
    %   minDecay=2;
    %   cmIx=floor(min(max((d(:,helpInd)-minDecay)/maxDecay,0),1)*127)+1;
    %   for dcIx=1:nDrugCondit
    %     curIx=drugID==ds.drugIx(dcIx);
    %     ph=scatter(d(curIx,pInd(1)),d(curIx,pInd(2)),3,cm(cmIx(curIx,:)),'filled');
    %     % set(ph,'MarkerEdgeColor',ds.pCol(dcIx,:));
    %   end
    %   set(gca,'xscale','log','yscale','log');
    %   axis([mima(:,1)' mima(:,2)'])
    %   set(gca,'xtick',anPSCPar{1,4},'ytick',anPSCPar{2,4})
    %   grid on
    %   xlabel(axLabel{1});
    %   ylabel(axLabel{2});
    %   box off
    %
    
