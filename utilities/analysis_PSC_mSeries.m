function [dCell1way,gCell1way,dCell2way,gCell2way]=analysis_PSC_mSeries(dataPath,plotPath,dataSet,fullPSCPar,ap)
% ** function [dCell1way,gCell1way,dCell2way,gCell2way]=analysis_PSC_mSeries(dataPath,dataSet,fullPSCPar,ap)
% generates a boxplot of normalized PSC parameter values for several data
% series (each processed by pscdeal.m) and puts out the numbers rearranged
% such that statistics on the differences between the series can easily be
% performed. The number of data series (e.g. different drugs) is not
% limited, but currently only one drug condition will be compared vs one
% other.
%
%                         >>> INPUT VARIABLES >>>
%
% NAME           TYPE/DEFAULT          DESCRIPTION
% dataPath       char                  path to data file
% plotPath       char                  path for figure files produced here
% dataSet        cell arr              list of data files etc. to be
%                                      processed (see template file)
% fullPSCPar     cell arr              list of parameters to be processed
%                                      (see template file)
% ap             struct                analysis/plot settings (see template
%                                      file)
%                     
%                         <<< OUTPUT VARIABLES <<<
%
% NAME           TYPE/DEFAULT           DESCRIPTION
%
%

% --- work in progress ---

nDs=size(dataSet,1);
nPar=size(fullPSCPar,1);

% cell containers for data for statistical analyses
dCell1way=cell(nPar,1);
dCell2way=cell(nPar,1);
% cell containers for group tag
gCell1way=cell(nPar,1);
gCell2way=cell(nPar,1);
% cell container for parameter and group tag (may be needed for boxplot)
gCellBoxplot=cell(nPar,1);


% ------ load & organize data ------
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
  
  for pIx=1:nPar
    
    parIx=strcmp(fullPSCPar{pIx,1},r.depPar);
    d=abs(r.pscrMn(:,indepParIx,parIx));
    
    % - kick all with any nan in first two columns
    d=d(all(isfinite(d(:,[1 2])),2),:);
    
    %     % ** simple pairwise nonparametric comparisons:
    %     stats=mes(d(:,1),d(:,2),'auroc','nBoot',1000);
    %     p=ranksum(d(:,1),d(:,2));
    %     disp([fullPSCPar{pIx,1} ': ' num2str(stats.auroc) ' [' num2str(stats.aurocCi') ']; p=' num2str(p)])
    
    % ** collecting data for mes2way: 
    n1=size(d,1);
    % - first and second column
    dCell2way{pIx}=cat(1,dCell2way{pIx},reshape(d(:,[1 2]),[n1*2 1]));
    % - group tag (factors: experimental series, drug condition)
    tmpGroup=[repmat(g,[n1*2 1]), reshape(ones(n1,1)*(1:2),[n1*2 1])];
    gCell2way{pIx}=cat(1,gCell2way{pIx},tmpGroup);
    
    % ** collecting data for mes1way and boxplot: normalize and keep only second column
    dCell1way{pIx}=cat(1,dCell1way{pIx},d(:,2)./d(:,1));
    % - group tag for mes1way (factor: experimental series)
    gCell1way{pIx}=cat(1,gCell1way{pIx},repmat(g,[n1 1]));
    % - group tag for boxplot (factors: parameter, experimental series)
    tmpGroup=[repmat(pIx,[n1,1]),repmat(g,[n1,1])];
    gCellBoxplot{pIx}=cat(1,gCellBoxplot{pIx},tmpGroup);
  end
end

% ------ boxplot ------
subplot(ap.subPlotGrid(1),ap.subPlotGrid(2),1)
ph=boxplot(cat(1,dCell1way{:}),cat(1,gCellBoxplot{:}),'widths',.9,'FactorGap',ap.factorGap,...
  'FactorSeparator',[1],'datalim',ap.dataLim,'colors',cat(1,dataSet{:,4}),...
  'whisker',ap.whisker,'symbol','+');
% undocumented Matlab: set properties of elements of boxplot
set(ph,'linewidth',1.0);
set(ph(1:2,:),'linestyle','-');
set(gca,'XTickLabel',{' '})

nicexy0ax(20);
% zero effect line
xl=get(gca,'xlim');
lh=line(xl,[1 1],'color','k','linestyle',':');
ylabel('norm. drug effect')
box off

% query whiskers' xData for proper placement of text
x=get(ph(1,:),'xdata');
x=cat(1,x{:});
x=x(:,1);
% compute mean of x position of boxes in each gruop
tmpArr=ones(nDs,1)*(1:nPar);
x=accumarray(tmpArr(:),x,[],@mean);
y=repmat(ap.dataLim(1)+diff(ap.dataLim)*.9,[nPar,1]);
th=text(x,y,fullPSCPar(:,2),'HorizontalAlignment','center');

% legend
subplot(ap.subPlotGrid(1),ap.subPlotGrid(2),2)
axis off
legend(gca,ph(1,1:nDs),dataSet(:,3));

% print?
if ~isempty(ap.printas)
  print(ap.printas,'-r300',[plotPath ap.plotFn])
end

