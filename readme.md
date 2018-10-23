# PSCFit

is Matlab code for the analysis of phasic and tonic postsynaptic currents (PSCs). Its core is a graphical user interface (GUI) for fitting exponentials to the decay time of phasic PSCs (pscfitgui). There is a number of additional primary and post-processing tools; see list below. The code makes use of in-house formats for time information ('time stamp lists', tsl).

Please note that the code in this repository is not self-sufficient, you'll additionally need the following repositories:
* fileIO
* etslfunc
* graphics
* sampledSeries
* utilities
 
Matlab toolboxes required:
* Signal Processing
* Statistics and Machine Learning
* Curve Fitting

### pscfitgui
GUI for fitting exponentials to the decay phase of PSCs. 

![snapshot](/doc/snapshot_pscfitgui.png)
Features:
* reading of raw data files (current traces) in Axon Binary Format (abf)
* reading of PSC time stamps (ideally produced by threshdetgui) from *.mat files
* optional pre-processing of the time series (lowpass and Savitzky-Golay filtering; up- and downsampling)
* flexible intervals for baseline computation, PSC peak search, and decay phase fitting 
* precise computation of rise time via upsampling of rising edge to 100 kHz
* fitting of one or two exponential components to the decay phase; computation of individual and weighted parameters (amplitude, decay time)
* flexible post-peak interval for fitting
* correction for stacked PSCs
* option to average cutouts and compute fit to average
* 2D scatter plots of parameters and plots of raw excerpts for inspection of events
* computation of three metrics of fit quality
* metric-based automatic and manual removal of individual events
* saving of all computed PSC parameters and analysis hyperparameters to file

#### utilities/pscdeal
Function which performs batch analysis and summary of PSC parameters. It collects parameters computed by pscfitgui, and/or computes a range of additional parameters, some of them based on detected PSCs, others based solely on the raw data. Specifically, the analyses and parameters are the following:
1. raw data-based: computation of the average base line level, base line noise, and time-averaged current carried by phasic PSCs; see function phantosic below
2. 'detected PSCs'-based: computation of parameters that can be extracted from the time stamps of the PSCs and the raw data, like rise time, amplitude and frequency of occurrence of PSCs
3. 'fitted PSCs'-based: collection and post-processing of parameters computed from the exponential fits to individual PSCs in pscfitgui, like the decay time

The function runs through a user-defined set of experiments, computes/collects all parameters, saves them in a file so that they can be easily accessed for plots and statistical analysis, and also optionally produces a summary figure. See template_call_pscdeal.m

![snapshot](/doc/snapshot_pscdeal.png)


#### utilities/phantosic
Function which computes the average base line level, base line noise, and the time-averaged current carried by phasic PSCs. Inspired by Glykys and Mody (J.Physiol. 582, 1163–1178, 2007). The function is called from within pscdeal, but can also be used independently.

![snapshot](/doc/snapshot_phantosic.png)

#### utilities/sumPlot_PSC
Function which generates collections of simple line plots of PSC parameters as processed by pscdeal. See template_call_sumPlot_PSC.m.

![snapshot](/doc/snapshot_sumPlot_PSC.png)


#### utilities/analysis_PSC_mSeries
Function useful for depicting and comparing PSC parameters across experimental series, each of them processed by pscdeal.
* generates a boxplot of normalized PSC parameter values
* puts out the numbers rearranged such that statistics on the differences between the series can easily be performed

See template_call_analysis_PSC_mSeries.m.

![snapshot](/doc/snapshot_analysis_PSC_mSeries.png)

#### utilities/analysis_IPSC_mvar
Performs multi-variate analysis of PSCs; experimental, under development. See template_call_analysis_IPSC_mvar.m.

![snapshot](/doc/snapshot_analysis_IPSC_mvar.png)


## General note on repositories in the ExpAnesth organization
The code in these repositories provides basic tools for the analysis of electrophysiological time series to members of the Section of Experimental Anesthesiology, Department of Anesthesiology, University Hospital of Tuebingen. Except where noted, code was written by Harald Hentschke. It has been designed primarily for in-house use by individuals who were instructed on its scope and limitations. Also, a substantial proportion of the code has been developed and extended over a time span of >10 years. In detail,

* the implementation of algorithms reflects the evolution of Matlab itself, that is, code that had been initially developed on older versions of Matlab does not necessarily feature newer techniques such as the new automatic array expansion as introduced in Matlab Release 2016b
* nonetheless, all code has been tested to run on Matlab R2018b
* while most m-files contain ample comments, documentation exists only for a few repositories
* checks of user input are implemented to varying degrees
* the code will be improved, updated and documented when and where the need arises