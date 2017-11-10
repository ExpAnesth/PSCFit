# Version history of PSCFitGui
## V1.0
First version of PSCFitGui
## V2.0
- cutouts of PSCs generated from abf file, no need anymore for the voluminous *cutout*.mat files produced by threshdetgui  
- PSCs can be deleted from within the scatter plot  
- option to specify sampling frequency and cutout window  
- click on main cutouts window toggles between continuous line plots and point-wise plots of cutouts  
## V 2.1 
- introduced two measures of fit error which can be accessed via the drop-down menus of the scatter plot
- the fitting algorithm now takes into account the recent history of the PSCS trace, enabling fits to PSCs riding on the decay phase of other PSCs; more generally speaking, to PSCs closely spaced in time  
- saving results automatically gets rid of events wich could not be fitted or even failed the preprocessing stage (previously this would only have been accomplished by hitting the 'delete events' button)  
- several minor bug/annoyance fixes  
## V 2.2
- new output parameter 'width' has been added (=half width of PSCs, as a complement or alternative to decay time constants)  
- PSCFit by default dumps both raw data and fits in the main workspace for easy post-analysis inspection of the results. See internal parameter wp.doDumpPSC.  
- renamed a number of core internal variables for consistency reasons  
- added option to lowpass filter the data prior to all analysis  
- bug fix: if a PSC's half width was more than the cutout length, an indexing error occurred, which was caught (including the misleading message display that this event could not be fitted). This resulted in slow events being completely omitted from analysis, whereas they should only have resulted in a missing value for half width. This is now corrected.  
- internal limit of slow decay component (ap.maxSlowDecay) set from 80 ms to 500 ms  
- extent of PSC to be fitted: now, a 'safety distance' to the following PSC is kept, thus ensuring that the following event's rise time is not regarded as part of the PSC  
- strings in some buttons colored  

## current version (V 2.3)
- improved PSC peak detection: now, a peak is looked for in window ranging from (approximately) xIntvPeak(1) to the minimum of [xIntvPeak(2), time to next PSC]; this prevents peaks from being missed when two PSCs follow each other closely  
- for more precise determination of rise time, rising flanks of IPSCs are now upsampled to 100 kHz  
- set wp.noiseHiCFreq to 1000 Hz (previously 500 Hz)  
- internal limit of slow decay component (ap.maxSlowDecay) set from 500 ms to 100 ms  
- numerous small fixes of code following code analyzer report  
- added utility functions, among them analysis_PSC_mSeries.m
- /utilities/pscdeal.m seriously renovated: i) global variables are gone, instead data are stored in userdata of a control GUI; ii) introduced new parameters, iii) reorganized parameters to be analyzed into three categories (from fitted PSCs, from raw data and tsl, from raw data only) and restructured and added code accordingly