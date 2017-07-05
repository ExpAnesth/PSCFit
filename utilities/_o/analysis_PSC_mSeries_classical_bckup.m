function analysis_PSC_mSeries_classical(dataPath,dataSet,fullPSCPar,ap)

nDs=size(dataSet,1);
nPar=size(fullPSCPar,1);

% cell container for data
dCell=cell(nPar,nDs);
% cell container for group tag
gCell=cell(nPar,nDs);

for g=1:nDs
  load([dataPath dataSet{g,1}], 'PSCRMN','ds','depPar');
  indepParIx=dataSet{g,2};
  
  for pIx=1:nPar
    parIx=find(strcmp(fullPSCPar{pIx,1},depPar));
    d=abs(PSCRMN(:,indepParIx,parIx));
    % normalize y first column
    d=d./repmat(d(:,1),[1 size(d,2)]);
    % ...not forgetting tag indicating parameter
    dCell(pIx,g)={d(:,2)};
    gCell(pIx,g)={repmat(pIx,[size(d,1) 1])};
  end
end


% plot 
deltaPar=2.5;
deltaDs=.7;

% define colors according to master file

mkfig(1,'b');
clf
subplot(2,2,1)
hold on

% loop over data sets 
ph=gobjects(1,nDs);
% § note: boxplot offers input args that allow specifying gaps between
% groups etc., so using a loop is unnecessary
for k=1:nDs
  xAx=(1:deltaPar:nPar*deltaPar)+(k-1)*deltaDs;
  tmpPh=boxplot(cat(1,dCell{:,k}),cat(1,gCell{:,k}),'positions',xAx,'widths',.5,...
    'datalim',[0 3],'colors',dataSet{k,4},'symbol','+');
  set(tmpPh,'linewidth',1.2);
  % don't ask
  ph(k)=tmpPh(6,1);
  % don't ask, either
  set(tmpPh(1:2,:),'linestyle','-')
end
nicexy0ax(15);
% zero effect line
xl=get(gca,'xlim');
lh=line(xl,[1 1],'color','k','linestyle',':');

set(gca,'xtick',(1:deltaPar:nPar*deltaPar)+(mean(1:nDs)-1)*deltaDs,...
  'xticklabel',fullPSCPar(:,2),'xticklabelrot',30);
ylabel('norm. drug effect')
box off

sph=subplot(2,2,2);
axis off
legend(sph,ph,dataSet(:,3));
