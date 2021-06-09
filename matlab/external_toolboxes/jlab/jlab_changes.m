%JLAB_CHANGES   Changes to JLAB in each release.
%
%   Changes new in version 1.5 
%
%   This is a major new release, which includes the following changes:
%       
%     --- Compatible with the graphics changes introduced in Matlab R2014b 
%     --- Completely new organization with a more modular approach
%     --- Many obsolete functions removed and redundancies eliminated
%
%   New functions
%
%   divgeom        - Geometric decomposition of eddy vorticity flux divergence.
%   blurspec       - Returns the blurred and aliased spectrum given the autocovariance.
%   findpath       - Returns the full pathname of a directory on the Matlab search path.
%   patchcontourf  - Generate filled contours using patches, with specified colors.
%   cellchunk      - Converts cell array data into uniform length 'chunks'.
%   inertialchunk  - Converts cell array data into chunks based on the Coriolis period.
%   cellmean       - Mean value of each element a cell array or set of cell arrays.
%   uv2latlon      - Integrates horizontal velocity to give latitude and longitude.
%   sig2latlon     - Converts an oscillatory signal to lat/lon displacements.
%   fminsearchbnd  - fminsearch, but with bound constraints by transformation. [By J. D'Errico]
%   whichdir       - Returns directory name containing file in search path.
%   jlab_allhelp   - Displays the help comments for all JLAB modules.
%
%   Improved datasets
%
%   drifters.mat   - Modified to include data through June 2014.  Please
%      note that the earlier verison of drifters.mat was missing half the 
%      data on account of a read error.  The new version corrects this.
%
%   Minor changes and improvements
%
%   SPHERESORT speed improvements, bug fix, and support for parallel computing.
%   READTOPO now works with Smith and Sandwell v. 18.1
%   FINDFILES now optionally searches directories recurvisely. 
%   PACKBOTH, PACKROWS, and PACKCOLS are now combined into PACKFIG.
%   JLAB_MAKEFIGS is now one file that makes figures for recent published papers. 
%   MATERNSPEC bugfix for input parameter arrays of length greater than one.
%   FLOATREGION now rejects individual trajectory segements outside of region.
%   CELLPLOT fixed minor bug in which CELLPLOT did not plot the first point.
%   TWOSPECPLOT no longer opens a new figure window by default.
%   RIDGEWALK now outputs ridges in a cell array, with one ridge per cell.  
%
%   ELLIPSEPLOT, MORSEBOX, FONTSIZE, TOPOPLOT, NOCONTOURS, TOPOPLOT, 
%   DISCRETECOLORBAR, JLAB_MAKEFIGS, and sample figures modified to work 
%   with Matlab version R2014b.
%
%   ELLVEL, ELLROSSBY, ELLBAND, ELLPARAMS, ELLRAD, and RIDGELEN now all 
%   correctly return empty arrays, or empty cell arrays, given empty input.
%
%   MATERNSPEC, MATERNCOV, MATERNIMP, and MATERNOISE improvements.  All 
%   extended to accommodate non-unit sample rates.
%
%   Removed a large number of obsolete functions; consolidated many others.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002--2015 J.M. Lilly --- type 'help jlab_license' for details

help jlab_changes

%   These were help back from version 0.95.
%
%   gaussprofile - Wavelet transform profile of a Gaussian with a Gaussian wavelet.
%   eddycensus   - Census of coherent eddies from streamfunction snapshots.
%   curveflow  Transport into and out of a region bounded by a closed curve.
%   eddyslice    - Slice through an eddy due to a turning background flow.
%   eddycurve    - Trajectory of an eddy center in a turning background flow.
%   eddyguess    - Eddy properties in a current meter record from a wavelet transform.
%   eddyfit      - Best fit parameters for an eddy advected past a current meter mooring.
%   wavetransderiv  - Rapidly calculate the wavelet transform of a signal's derivative.