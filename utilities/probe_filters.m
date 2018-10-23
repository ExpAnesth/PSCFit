% This is a collection of ad-hoc code cells illustrating the effects of
% varying preprocessing filters on the signals.

%% test lowpass filters
loCFreq=[nan 3000 1000 15]';
nFreq=numel(loCFreq);
[d,si]=abfload('e:\_data\_IPSCFit\2015_10_16_0011_exc.abf');
% time axis in ms
t=(1:numel(d))*(si/1000);
figure(1)
clf
ph=gobjects(nFreq,1);
hold on
for ii=1:nFreq
  if isnan(loCFreq(ii))
    % unfiltered data
    ph(ii)=plot(t,d);
    set(ph(ii),'color',[.7 .7 .7])
  elseif isfinite(loCFreq(ii)) && loCFreq(ii)>100
    % butterworth
    dFi=lofi(d,si,loCFreq(ii));
    ph(ii)=plot(t,dFi);
  else
    % Savitzky-Golay filter
    dFi=sgolayfilt(d,5,loCFreq(ii));
    ph(ii)=plot(t,dFi);
  end
  
end
legend(ph,int2str(loCFreq))

%% test highpass filters (for separation of noise from signal
hiCFreq=[nan 500 1000]';
nFreq=numel(hiCFreq);
[d,si]=abfload('e:\_data\_IPSCFit\2015_10_16_0011_exc.abf');
% time axis in ms
t=(1:numel(d))*(si/1000);
figure(1)
clf
ph=gobjects(nFreq,1);
hold on
for ii=1:nFreq
  if isnan(hiCFreq(ii))
    % unfiltered data
    ph(ii)=plot(t,d+500);
    set(ph(ii),'color',[.7 .7 .7])
  elseif isfinite(hiCFreq(ii)) && hiCFreq(ii)>100
    % butterworth
    dFi=hifi(d,si,hiCFreq(ii));
    ph(ii)=plot(t,dFi);
%   else
%     % Savitzky-Golay filter
%     dFi=sgolayfilt(d,5,hiCFreq(ii));
%     ph(ii)=plot(t,dFi);
  end
  
end
legend(ph,int2str(hiCFreq))
