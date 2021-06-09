% jMatern:  Parametric spectral analysis based on the Matern process
%
% Top-level functions
%   maternoise - Generates realizations of the Matern random process.  [with A. Sykulski]      
%   maternspec - Fourier spectrum of the Matern random process.                                
%   materncov  - Autocovariance function of the Matern random process.                         
%   maternimp  - Impulse response function for the Matern random process.                      
%
% Other utilities
%   blurspec   - Returns the blurred and aliased spectrum given the autocovariance.
%   fminsearchbnd: - FMINSEARCH, but with bound constraints by transformation. [By J. D'Errico]
%
% Low-level Matern functions
%   maternc    - Returns the normalization function C_ALPHA for a Matern process.
%   maternchol - Cholesky decomposition of Matern the covariance matrix.
%   maternedge - Long-time cutoff edge for the Matern impulse response function.               
%   maternfun  - Returns the Matern function.                                                  
%
% See also jSpectral.

help jMatern

if 0
% Top-level functions
   maternoise %- Generates realizations of the Matern random process.  [with A. Sykulski]      
   maternspec %- Fourier spectrum of the Matern random process.   
   materncov  %- Autocovariance function of the Matern random process.                         
   maternimp  %- Impulse response function for the Matern random process.                      

   % Other utilities
   blurspec   %- Returns the blurred and aliased spectrum given the autocovariance.
   fminsearchbnd: %- FMINSEARCH, but with bound constraints by transformation. [By J. D'Errico]

%   Low-level functions
   maternc    %- Returns the normalization function C_ALPHA for a Matern process.
   maternchol %- Cholesky decomposition of Matern covariance matrix.
   maternedge %- Long-time cutoff edge for the Matern impulse response function.               
   maternfun  %- Returns the Matern function.                                                  
end
