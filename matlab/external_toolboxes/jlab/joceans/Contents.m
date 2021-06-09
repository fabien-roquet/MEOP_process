% jOceans:  Oceanography-specific data and model analysis tools
%
% Conversions for Lagrangian trajectories
%   latlon2uv   - Converts latitude and longitude to horizontal velocity.  
%   uv2latlon   - Integrates horizontal velocity to give latitude and longitude.  
%
% Manipulating Lagrangian trajectories
%   trajextract - Extracts Lagrangian trajectory segments within given region.        
%   trajfill    - Fills float or drifter trajectories with linear interpolation.       
%   trajunwrap  - Unwraps Lagrangian trajectories from a periodic domain.              
%   trajwrap    - Wraps Lagrangian trajectories to fit within a periodic domain.       
%   trajchunk   - Converts cell array data into chunks based on the Coriolis period.
%
% Idealized numerical model tools
%   psi2fields   - Velocity and other fields from the streamfunction. [with P.E. Isachsen]    
%   periodize    - Returns a doubly periodic version of an input array.   
%
% Eulerian eddy identification and analysis
%   inellipse    - Locates points on the interior of ellipses.                          
%   closedcurves - Locate and interpolate closed curves in a possibly periodic domain.
%   curvemoments - Centroid, area, and many other moments of a closed curve.          
%   divgeom      - Geometric decomposition of eddy vorticity flux divergence.     
%
%  Plotting tools for mooring data
%   hodograph  - Generate hodograph plots (simple and fancy).                                   
%   provec     - Generate progressive vector diagrams (simple and fancy).             
%   stickvect  - Plots "stick vectors" for multicomponent velocity time series.       
%
%  See also jSphere, jMatern.

%   Low-level functions
%   curveinterp  - Interpolate a field or its gradient onto a set of curves.           
%   orbitbreaks   - Separate orbit into passes based on turning points.       

help joceans
if 0
%Conversions for Lagrangian trajectories
  latlon2uv   %- Converts latitude and longitude to horizontal velocity.  
  uv2latlon   %- Integrates horizontal velocity to give latitude and longitude.  
%  [See also jSphere] 

%Manipulating Lagrangian trajectories
  trajextract %- Extracts Lagrangian trajectory segments within given region.        
  trajfill    %- Fills float or drifter trajectories with linear interpolation.       
  trajunwrap  %- Unwraps Lagrangian trajectories from a periodic domain.              
  trajwrap    %- Wraps Lagrangian trajectories to fit within a periodic domain.       
  trajchunk   %- Converts cell array data into chunks based on the Coriolis period.

%Idealized numerical model tools
  psi2fields   %- Velocity, vorticity, and strain fields from the streamfunction.      
  periodize    %- Returns a doubly periodic version of an input array.   

%Eulerian eddy identification and analysis
  inellipse    %- Locates points on the interior of ellipses.                          
  closedcurves %- Locate and interpolate closed curves in a possibly periodic domain.
  curvemoments %- Centroid, area, and many other moments of a closed curve.          
  divgeom      %- Geometric decomposition of eddy vorticity flux divergence.     

% Plotting tools
  hodograph  %- Generate hodograph plots (simple and fancy).                                   
  provec     %- Generate progressive vector diagrams (simple and fancy).             
  stickvect  %- Plots "stick vectors" for multicomponent velocity time series.       

%  Low-level functions
  curveinterp  %- Interpolate a field or its gradient onto a set of curves.           
  orbitbreaks  %- Separate orbit into passes based on turning points. 
end
