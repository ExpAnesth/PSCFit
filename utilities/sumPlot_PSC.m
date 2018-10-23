function sumPlot_PSC(dataPath,plotPath,dataSet,fullPSCPar,ap)
% ** function sumPlot_PSC(dataPath,plotPath,dataSet,fullPSCPar,ap)
% produces simple line plots of IPSC parameters, one figure per data set 

nDs=size(dataSet,1);
nPar=size(fullPSCPar,1);

for g=1:nDs
  % load results, differentiating between old (individual variables) and
  % new (struct r) styles
  v=whos('-file',[dataPath dataSet{g,1}]);
  if ~isempty(intersect({v.name},'r'))
    load([dataPath dataSet{g,1}], 'r');
  else
    load([dataPath dataSet{g,1}], 'PSCRMN','depPar');
    % make fields of r
    r.pscrMn=PSCRMN;
    r.depPar=depPar;
    clear PSCRMN depPar
  end
  indepParIx=dataSet{g,2};
  xlab=dataSet{g,3};
  
  mkfig(g,'b');
  clf
  for pIx=1:nPar
    % slice index into r.pscrMn
    parIx=strcmp(fullPSCPar{pIx,1},r.depPar);
    if ~any(parIx)
      error(['parameter ''' fullPSCPar{pIx,1} ''' not present in data set ' dataSet{g,1}])
    end
    
    % --- data retrieval
    d=abs(r.pscrMn(:,indepParIx,parIx));

    % --- plotting part
    subplot(ap.subPlotGrid(1),ap.subPlotGrid(2),pIx)
    ylab=fullPSCPar{pIx,2};
    yscl=fullPSCPar{pIx,3};
    [ph,~]=avplot(abs(r.pscrMn(:,indepParIx,parIx)),'x',1:numel(indepParIx),...
      'avType','md');
    set(gca,'yscale',yscl,'xtick',1:numel(indepParIx),'xticklabel',xlab,...
      'xticklabelrot',30);
    set(ph,'color',[.5 .5 .5]);
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

