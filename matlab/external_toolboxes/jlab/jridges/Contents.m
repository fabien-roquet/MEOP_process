% jRidges:  Wavelet ridge analysis of modulated oscillatory signals
%  
%  Top-level functions
%   ridgewalk  - Extract wavelet transform ridges, including bias estimates. 
%   ridgemap   - Maps ridge quantities back onto the time series.            
%   instmom    - Univariate and multivariate instantaneous moments.          
%
%  Ridge utilities
%   ridgelen   - Wavelet ridge length expressed as number of full cycles.    
%   periodindex - Returns time index in increments of instantaneous period.  
%
%  See also jEllipse, jWavelet.

%   Low-level functions
%   isridgepoint - Finds wavelet ridge points using one of several criterion.
%   lininterp  - Fast linear interpolation for arbitrary-sized arrays.       
%   quadinterp - Fast quadratic interpolation for arbitrary-sized arrays.    
%   ridgechains - Forms ridge curves by connecting transform ridge points.   
%   ridgeinterp - Interpolate quantity values onto ridge locations.          

help jRidges
if 0
%jRidges:  Wavelet ridge analysis
 
% Top-level functions
  ridgewalk  %- Extract wavelet transform ridges, including bias estimates. 
  ridgemap   %- Maps ridge quantities back onto the time series.            
  instmom    %- Univariate and multivariate instantaneous moments.          

% Ridge utilities
  ridgelen   %- Wavelet ridge length expressed as number of full cycles.    
  periodindex %- Returns time index in increments of instantaneous period.  

% See also jEllipse, jWavelet.

  isridgepoint %- Finds wavelet ridge points using one of several criterion.
  lininterp  %- Fast linear interpolation for arbitrary-sized arrays.       
  quadinterp %- Fast quadratic interpolation for arbitrary-sized arrays.    
  ridgechains %- Forms ridge curves by connecting transform ridge points.   
  ridgeinterp %- Interpolate quantity values onto ridge locations.          
end
