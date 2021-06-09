%  jSpectral:  Multitaper spectral analysis, and other time series tools
% 
%  Multitaper spectral analysis
%   sleptap    - Calculate Slepian tapers.                                        
%   mspec      - Multitaper power and cross spectra.       
%
%  Multitaper polarization analysis
%   msvd       - Singular value decomposition for polarization analysis.   
%   polparams  - Spectral matrix polarization parameters.                         
%   specdiag   - Diagonalize a 2 x 2 spectral matrix.        
%
% Assorted other transforms 
%   slidetrans  - Sliding-window ('moving-window') Fourier transform.   
%   anatrans    - Analytic part of signal.                                         
%   wigdist     - Wigner distribtion (alias-free algorithm).   
% 
% Time series analysis utilities
%   doublen     - Interpolates a time series to double its length.                 
%   fourier     - The one-sided Fourier frequencies for a given length time series.
%   sampletimes - Computes mean sampling intervals and their statistics.          
%
%  Plotting tools
%   twospecplot - Plots a pair of rotary or Cartesian spectra.   
%
%  See also jWavelet, jEllipse, jMatern.

%   Low-level functions
%   timeseries_boundary - Apply boundary conditions to data before transform.


help jspectral

if 0
% jSpectral:  Multitaper spectral analysis, and other time series tools

% Multitaper spectral analysis
  sleptap    %- Calculate Slepian tapers.                                        
  mspec      %- Multitaper power and cross spectra.       

% Multitaper polarization analysis
  msvd       %- Singular value decomposition for polarization analysis.   
  polparams  %- Spectral matrix polarization parameters.                         
  specdiag   %- Diagonalize a 2 x 2 spectral matrix.        

% Assorted other transforms 
  slidetrans  %- Sliding-window ('moving-window') Fourier transform.   
  anatrans    %- Analytic part of signal.                                         
  wigdist     %- Wigner distribtion (alias-free algorithm).   

% Time series analysis utilities
  doublen     %- Interpolates a time series to double its length.                 
  fourier     %- The one-sided Fourier frequencies for a given length time series.
  sampletimes %- Computes mean sampling intervals and their statistics.          

% Plotting tools
  twospecplot %- Plots a pair of rotary or Cartesian spectra.   
  
% Low-level functions
  timeseries_boundary %- Apply boundary conditions to data before transform.

% See also jWavelet, jEllipse, jMatern.      
end
