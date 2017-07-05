function sumPlot_PSC(dataPath,plotPath,dataSet,fullPSCPar,ap)
% produces simple line plots of IPSC parameters 


nDs=size(dataSet,1);
nPar=size(fullPSCPar,1);


for g=1:nDs
  % plot
  mkfig(g,'b');
  clf
  
  load([dataPath dataSet{g,1}], 'PSCRMN','ds','depPar','LISTEXP');
  indepParIx=dataSet{g,2};
  xlab=dataSet{g,3};
  
  for pIx=1:nPar
    % slice index into PSCRMN
    parIx=find(strcmp(fullPSCPar{pIx,1},depPar));
    
    % --- data retrieval
    d=abs(PSCRMN(:,indepParIx,parIx));
    %     % - kick all with any nan in first two columns
    %     tmpD=d(all(isfinite(d(:,[1 2])),2),:);
    %
    %     % ** simple pairwise nonparametric comparisons:
    %     stats=mes(tmpD(:,1),tmpD(:,2),'auroc','nBoot',1000);
    %     p=ranksum(tmpD(:,1),tmpD(:,2));
    %     disp([fullPSCPar{pIx,1} ': ' num2str(stats.auroc) ' [' num2str(stats.aurocCi') ']; p=' num2str(p)])

    % --- plotting part
    subplot(ap.subPlotGrid(1),ap.subPlotGrid(2),pIx)
    ylab=fullPSCPar{pIx,2};
    yscl=fullPSCPar{pIx,3};
    [ph,ebh]=avplot(abs(PSCRMN(:,indepParIx,parIx)),'x',1:numel(indepParIx),...
      'avType','md');
    set(gca,'yscale',yscl,'xtick',ds.indepPar(indepParIx),'xticklabel',xlab,...
      'xticklabelrot',30);
    set(ph,'color',[.5 .5 .5]);
    % set(ebh,'markersize',13,'markerfacecolor','k','linewidth',3);
    grid on
    ylabel(ylab)
    nicexyax(12);
    
  end
  drawnow
  
  % print?
  if ~isempty(ap.printas)
    print(ap.printas,'-r300',[plotPath strrep(dataSet{g,1},'.mat','_paramSum')])
  end
  
end

