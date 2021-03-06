%JLAB_HIGHLIGHTS  Introduction to some of the most useful routines.
%
%   JLAB_HIGHLIGHTS  Highlights from the JLAB toolbox.
%
%   TWODSTATS, TWODMED, TWODHIST.  Examining temperature mean and standard
%   deviation over a lat/lon grid, or current speed as a function of say
%   temperature and salinity, or any analysis of multivariate datasets, is
%   faster and easier than you ever imagined.  If you are working with 
%   data and not using these routines, you should seriously take five 
%   minutes right now and try them out.  I wrote fast, exact, and loopless
%   algorithms that are often 100 times (!) faster than the obvious way.  
%
%   MSPEC.  If you analyze time series, you need a method you can rely on.
%   There are many compelling reasons to prefer the "multitaper" method.
%   In this method, a set of estimates from optimally concentrated tapers
%   or window functions are averaged together, reducing variance while 
%   minimizing bias.  MSPEC makes using the multitaper straightforward.  
%   Just call SLEPTAP first to generate the tapers.  MSPEC features support
%   for cross spectra and rotary spectra, as well as 'adaptive' estimation.
%   
%   WAVETRANS.  A lot of work and original research has gone into this 
%   continuous wavelet transform routine, which I believe is the best you 
%   will find anywhere, and it's free.  Features natural integration with
%   the generalized Morse wavelets (see below), choice of conditions at 
%   the signal endpoints, support for multiple wavelets or mulit-component
%   data, and convenient handling of real-valued or complex-valued data.
%
%   MORSEWAVE.  A detailed analysis of types of continuous wavelets 
%   (Lilly and Olhede, 2009) points to a super-family, the generalized 
%   Morse wavelets, that should answer your questions about which wavelet
%   to employ. Pick the generalized Morse wavelet with parameter gamma=3 
%   unless you have a very good reason to do otherwise.  This wavelet is
%   'like' the Morlet in spirit, but without the Morlet's serious flaws at
%   narrow time-domain settings---which is when you most want to employ a 
%   wavelet anyway.  Use MORSESPACE to easily determine frequency bins.
%
%   RIDGEWALK.  If your time series data contains signals that oscillate
%   but also change in time, the wavelet transform can be used as the basis
%   for a remarkably powerful analysis that connects time-varying 
%   properties to spectral structure.  The starting place for this analysis
%   is RIDGEWALK.  This advanced algorithm incorporates original research 
%   from Lilly & Olhede (2010a+b, 2011).  A key innovation is the ability 
%   to view *multivariate* and univariate time series within the same 
%   framework as modulated oscillations.
%
%   POLYSMOOTH.  Need to map scattered data, for example, temperature as a 
%   function of X and Y?  POLYSMOOTH, short for "Polynomial Smoothing", 
%   is a powerful and quite general mapping algorithm.  Features choice of 
%   order of the fit, choice of local weighting function, support for 
%   spherical geometry, options for constant fit radius or constant number
%   of data points with variable radius.  This algorithm was implemented
%   for the forthcoming Aquarius Sea Surface Salinity satellite.
%
%   JGRAPH.  Finessing figures in Matlab can be a total hassle.  Which
%   would you prefer to type, version (a) or  version (b) ? 
%
%      (a)  z=peaks;
%           plot(z(:,1),'linewidth',3,'linestyle','-','color','b'),hold on
%           plot(z(:,2),'linewidth',2,'linestyle','--','color','r')
%
%      (b)  z=peaks;plot(z(:,1:2)), linestyle 3b- 2r--
%
%   They are identical.  JLAB's LINESTYLE will save you a lot of time.  
%   It includes natural syntax for grayscale, as in "linestyle 3k- 2G-.",
%   where the letter G stands for the shade of gray, between A and K. 
%   The JGRAPH module contains many similar time-saving functions, such as
%   YLOG and YLIN, YOFFSET, FLIPY, OUTTICKS, FLIPMAP, etc.
%
%   PACKFIG.  Are your multiple subplots too far away from each other?  
%   Don't like having to set the axes properties by hand each time?  Try  
%   this function.  For example, "packfig(3,4)" will squish your 3 rows and 
%   4 columns so they are right next to each other, removing the redundant
%   x- and y-axes labels.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2015 J.M. Lilly --- type 'help jlab_license' for details

