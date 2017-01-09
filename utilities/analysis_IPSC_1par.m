% a collection of various plots of IPSC parameters  - pscdeal must have run
% before and global variables PSCRMN and PSCR containing the values as well
% as struct ds and variable depPar must exist

%% plot of single parameter vs drug application
labelscale('fontSz',16,'scaleFac',.5,'lineW',2,'markSz',8); 
figure(10), 
clf
hold on

parIx=find(strcmp('tDecay',depPar));
ylab='{\tau}_{decay} (ms)';

parIx=find(strcmp('chargePscTot',depPar));
ylab='charge/time (pA)';

parIx=find(strcmp('freq',depPar));
ylab='frequency (Hz)';

% indexes into PSCRMN, corresponding to drug conditions listed in
% ds.indepPar and ds.indepParLabel
indepParIx=[1 2 3];

% % define colors corresponding to GABAAR profile
% fac=graphicsDef_ACh;
% fIx=find(strcmp('GABAAR profile',{fac.name}));
% if 1
%   % zol
%   colIx=find(strcmp('a1',fac(fIx).levelName));
% else
%   % dia in WT
%   colIx=find(strcmp('a1235',fac(fIx).levelName));
% end
% % finally, the color
% col=fac(fIx).color(colIx,:);
col='k';

[ph,ebh]=avplot(PSCRMN(:,indepParIx,parIx),'x',ds.indepPar(indepParIx),'mType','v');
set(gca,'xtick',ds.indepPar(indepParIx),...
  'xticklabel',ds.indepParLabel(indepParIx));
set(ph,'linewidth',1,'color','k');
set(ebh,'markersize',24,'markerfacecolor','k','linewidth',3);
ylabel(ylab)
nicexyax(10);

%% plot of one par vs other (mean values)
% §§ make stats & their plots nonparametric
labelscale('fontSz',16,'scaleFac',.5,'lineW',2,'markSz',8); 
figure(11), clf, hold on

parIx=[find(strcmp('freq',depPar)) find(strcmp('amp',depPar))];
xlab='frequency (Hz)';
ylab='peak amplitude (pA)';

% indexes into PSCRMN, corresponding to drug conditions listed in
% ds.indepPar and ds.indepParLabel
indepParIx=[2 3];

% define colors corresponding to ACh status 
fac=graphicsDef_ACh;
fIx=find(strcmp('ACh status',{fac.name}));
col1Ix=find(strcmp('ACh+',fac(fIx).levelName));
col2Ix=find(strcmp('ACh-',fac(fIx).levelName));
% finally, the colors
col1=fac(fIx).color(col1Ix,:);
col2=fac(fIx).color(col2Ix,:);

x1=PSCRMN(:,indepParIx(1),parIx(1));
y1=PSCRMN(:,indepParIx(1),parIx(2));
x2=PSCRMN(:,indepParIx(2),parIx(1));
y2=PSCRMN(:,indepParIx(2),parIx(2));

% line
ph=plot(PSCRMN(:,indepParIx,parIx(1))',PSCRMN(:,indepParIx,parIx(2))','-');
set(ph,'color',[0 0 0],'linewidth',1);
% condition 1
ph=plot(x1,y1,'kv');
set(ph,'markerfacecolor','w','color',col1);
% condition 2
ph=plot(x2,y2,'kv')
set(ph,'markerfacecolor','w','color',col2);

% condition 1: mean & std and plot of average with 95% CI
mdx1=mean(x1);
stx1=std(x1);
nx1=numel(x1);
semx1=stx1./sqrt(nx1);
cix1=semx1*tinv(.05/2,nx1);

mdy1=mean(y1);
sty1=std(y1);
ny1=numel(y1);
semy1=sty1./sqrt(ny1);
ciy1=semy1*tinv(.05/2,ny1);

lh=errorcross2(...
  [mdx1,mdy1],...
  [cix1 ciy1],...
  [cix1 ciy1],...
  'color',col1,'linewidth',3);

ph=plot(mdx1,mdy1,'kv');
set(ph,'markersize',24,'markerfacecolor',col1,'linewidth',3);


% same procedure for condition2
mdx2=mean(x2);
stx2=std(x2);
nx2=numel(x2);
semx2=stx2./sqrt(nx2);
cix2=semx2*tinv(.05/2,nx2);

mdy2=mean(y2);
sty2=std(y2);
ny2=numel(y2);
semy2=sty2./sqrt(ny2);
ciy2=semy2*tinv(.05/2,ny2);

lh=errorcross2(...
  [mdx2,mdy2],...
  [cix2 ciy2],...
  [cix2 ciy2],...
  'color',col2,'linewidth',3);

ph=plot(mdx2,mdy2,'kv');
set(ph,'markersize',24,'markerfacecolor',col2,'linewidth',3);

xlabel(xlab);
ylabel(ylab);

nicexyax;

%% joint IEI of IPSCs

%% comparison of amp vs ampRise (one conc)

figure(1), 
clf
hold on

anPar={'allAmp','amp'};
ylab='peak amplitude (pA)';

% indexes into PSCR, corresponding to drug conditions listed in
% ds.indepPar and ds.indepParLabel
indepParIx=[2];

col='k';

for g=1:numel(LISTEXP)
  d1=PSCR{g,indepParIx,strcmp(anPar{1},depPar)};
  d2=PSCR{g,indepParIx,strcmp(anPar{2},depPar)};
  [n1,edges]=histcounts(d1,40);
  n2=histcounts(d2,edges);
  subplot(3,4,g)
  hold on
  % ph=stairs(edges(1:end-1),[n1' n2']);
  % plot ratio
  ph=stairs(edges(1:end-1),[n2'./ n1']);
  axis tight
  grid on
  
end

%% comparison of fitted vs raw parameters, e.g. chargePscTot with chargePhas
anPar={'chargePscTot','chargePhas'};
% anPar={'freqFit','freq'};
anPar={'allAmp','chargePhas'};


figure(1), 
clf
subplot(3,1,1)
plot(PSCRMN(:,:,strcmp(anPar{1},depPar)),PSCRMN(:,:,strcmp(anPar{2},depPar)),'o');
nicexyax;
grid on
xlabel(anPar{1});
ylabel(anPar{2});

subplot(3,1,2)
avplot(PSCRMN(:,:,strcmp(anPar{1},depPar))./PSCRMN(:,:,strcmp(anPar{2},depPar)));
xlabel([anPar{1} '/' anPar{2}])


