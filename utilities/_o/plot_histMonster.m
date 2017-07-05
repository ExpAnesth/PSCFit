%%%PSCdeal needs to run first and wite variable HISTMONSTER. Here, this
%%%variable is averaged in all 3 dimensions (amp x tdecay x frequency) and
%%%one cumulative plot is created. Frequency is normalized to the first
%%%drug condition (or the condition specified). Also, highest peak across
%%%all drug conditions is being determined and colors are adjusted
%%%according to highest peak.


figure(1), clf
nSp=4;
nCol=2;
nRow=ceil(nSp/nCol);

% index of element (in fourth dim) in HISTMONSTER representing control
% condition (set to Nan for no normalization)
normIx=1;
if isfinite(normIx)
  n=HISTMONSTER(:,:,:,normIx);
  normFac=repmat(sum(sum(n)),size(n,1),size(n,2));
end

cLim=[];
sph=[];
for g=1:size(HISTMONSTER,4)
  n=HISTMONSTER(:,:,:,g);
  % normalize (number of events)?
  if isfinite(normIx)
    n=n./normFac;
  end
  % kick out entries of all zeros
  n=n(:,:,any(any(n)));
  % average
  n=nanmedian(n,3);
  cLim(g)=prctile(n(:),98);
  sph(g)=subplot(nRow,nCol,g);
  [~,h]=contourf(n,30);
  set(h,'linestyle','none');
  grid on
  ylabel(ds.plotPar1); xlabel(ds.plotPar2);
  title(ds.indepParLabel{g});
  % §§ shitty
  xtick=get(gca,'xtick');
  ytick=get(gca,'ytick');
  set(gca,'xtick',xtick,'xticklabel',int2str(round(ds.plotPar2Bin(xtick)')));
  set(gca,'ytick',ytick,'yticklabel',int2str(round(ds.plotPar1Bin(ytick)')));
end
set(sph,'clim',[0 max(cLim)]);