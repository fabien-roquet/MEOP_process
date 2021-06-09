% jCommon:  Useful general-purpose functions, common to other JLAB modules
%
% Array functions
%   aresame    - Test whether two N-D arrays are the same.                        
%   arrayify   - Converts a set of scalars or arrays into column arrays.            
%   blocklen   - Counts the lengths of 'blocks' in an array.                        
%   blocknum   - Numbers the contiguous blocks of an array.   
%   lnsd       - Last non-singleton dimension of an array.                          
%   matmult    - Matrix multiplication for arrays of matrices.                      
%   nonnan     - Return all non-NAN elements of an array.                           
%   vectmult   - Matrix multiplication for arrays of vectors.    
%
% File and directory tools
%   commentlines - Returns the comment lines from m-files.                          
%   findfiles    - Returns all files in a directory with a specified extension.       
%   findpath     - Returns the full pathname of a directory on the Matlab search path.
%   jhelp        - Opens linked JLAB help files in Matlab's internal web browser.     
%   whichdir     - Returns directory name containing file in search path.             
%
% Mathematical aliases
%   choose     - Binomial coefficient: CHOOSE(N,K) = "N choose K"
%   frac       - Fraction: FRAC(A,B)=A./B                                                    
%   imlog      - Imaginary part of log: IMLOG(X)=UNWRAP(IMAG(LOG(X)))                                      
%   iseven     - True for even integer values; false otherwise.                     
%   isodd      - True for odd integer values; false otherwise. 
%   oprod      - Outer product:  OPROD(X,Y)=X*CONJ(Y')                              
%   res        - Residual after flooring:  RES(X)=X-FLOOR(X)                        
%   rot        - Complex-valued rotation:  ROT(X)=EXP(SQRT(-1)*X)  
%   squared    - Squares the modulus of its argument:  SQUARED(X)=ABS(X).^2         
% 
% Matrices, polynomials, and special functions
%   jmat       - 2x2 rotation matrix through specified angle.                       
%   jmat3      - 3x3 rotation matrix through specified angle.  
%   tmat       - 2x2 complex grouping matrix.  TMAT = [1  i; 1 -i] / SQRT(2)
%   bellpoly   - Complete Bell polynomials.  
%   hermpoly   - Hermite polynomials. [with F. Rekibi]   
%   hermfun    - Orthonormal Hermite functions. [with F. Rekibi]                    
%
% Dataset organization as structures
%   make       - Create a structure containing named variables as fields.  
%   matsave    - Create and save structure of variables as a mat-file. 
%   use        - Copies structure fields into named variables in workspace.         
%
% Statistics
%   cum2mom    - Convert cumulants to moments.    
%   mom2cum    - Convert moments to cumulants.                                      
%   pdfprops   - Mean and variance associated with a probability distribution.  
%   simplepdf  - Gaussian, uniform, Cauchy, and exponential pdfs.                   
%
% Filling bad data points
%   fillbad    - Linearly interpolate over bad data points.  
           
%   Low-level functions
%   reporttest - Reports the result of an m-file function auto-test.                
%   to_grab_from_caller - Returns a string to grab variable values from caller.     
%   to_overwrite        - Returns a string to overwrite original arguments.  

help jCommon
if 0
% jCommon:  General-purpose functions, used by other JLAB modules
%
% Array functions
   aresame    %- Test whether two N-D arrays are the same.                        
   arrayify   %- Converts a set of scalars or arrays into column arrays.            
   blocklen   %- Counts the lengths of 'blocks' in an array.                        
   blocknum   %- Numbers the contiguous blocks of an array.   
   lnsd       %- Last non-singleton dimension of an array.                          
   matmult    %- Matrix multiplication for arrays of matrices.                      
   nonnan     %- Return all non-NAN elements of an array.                           
   vectmult   %- Matrix multiplication for arrays of vectors.    
%
% File and directory tools
   commentlines %- Returns the comment lines from m-files.                          
   findfiles    %- Returns all files in a directory with a specified extension.       
   findpath     %- Returns the full pathname of a directory on the Matlab search path.
   jhelp        %- Opens linked JLAB help files in Matlab's internal web browser.     
   whichdir     %- Returns directory name containing file in search path.             
%
% Mathematical aliases
   choose     %- Binomial coefficient: CHOOSE(N,K) = "N choose K"  
   frac       %- FRAC(A,B)=A./B;                                                    
   imlog      %- IMLOG(X)=UNWRAP(IMAG(LOG(X)))                                      
   iseven     %- True for even integer values; false otherwise.                     
   isodd      %- True for odd integer values; false otherwise. 
   oprod      %- Outer product:  OPROD(X,Y)=X*CONJ(Y')                              
   res        %- Residual after flooring.  RES(X)=X-FLOOR(X)                        
   rot        %- Complex-valued rotation:  ROT(X)=EXP(SQRT(-1)*X)  
   squared    %- Squares the modulus of its argument:  SQUARED(X)=ABS(X).^2         
% 
% Matrices, polynomials, and special functions
   jmat       %- 2x2 rotation matrix through specified angle.                       
   jmat3      %- 3x3 rotation matrix through specified angle.  
   tmat       %- 2x2 complex grouping matrix.  TMAT = [1  i; 1 -i] / SQRT(2).  
   bellpoly   %- Complete Bell polynomials.  
   hermpoly   %- Hermite polynomials. [with F. Rekibi]   
   hermfun    %- Orthonormal Hermite functions. [with F. Rekibi]                    
%
% Dataset organization as structures
   make       %- Create a structure containing named variables as fields.  
   matsave    %- Create and save structure of variables as a mat-file. 
   use        %- Copies structure fields into named variables in workspace.         
%
% Statistics
   cum2mom    %- Convert cumulants to moments.    
   mom2cum    %- Convert moments to cumulants.                                      
   pdfprops   %- Mean and variance associated with a probability distribution.  
   simplepdf  %- Gaussian, uniform, Cauchy, and exponential pdfs.                   
%
% Filling bad data points
   fillbad    %- Linearly interpolate over bad data points.  
           
%   Low-level functions
   reporttest %- Reports the result of an m-file function auto-test.                
   to_grab_from_caller %- Returns a string to grab variable values from caller.     
   to_overwrite        %- Returns a string to overwrite original arguments.  
end