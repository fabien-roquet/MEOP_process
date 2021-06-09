function[varargout]=col2cell(varargin)
%COL2CELL  Converts 'column-appended' data into cell arrays of column vectors.
%
%   X=COL2CELL(COL) converts the column vector COL of N data blocks into a
%   cell array X containing N column vectors.
%
%   The column vector COL has a NAN marking the end of each block.  The
%   cell array X exludes the NANs.  
%
%   Since NANs mark the end of each block, missing or bad data in COL
%   should instead be indicated by INFs.
%   __________________________________________________________________
%
%   Multiple input /output arguments
%
%   [X1,X2,...,XN]=COL2CELL(C1,C2,...,CN) also works, where the CN are all 
%   the same size.
% 
%   In this case C1 is the 'key' and the locations of NANs in C1 is used 
%   to break all the CN into cells.   Thus the CN with N>1 may contain
%   additional NANS.
%   __________________________________________________________________
% 
%   Variable overwriting
%
%   COL2CELL(C1,C2,...,CN) with no output arguments overwrites the input
%   variables.
%   __________________________________________________________________
%
%   Invertibility
%
%   COL2CELL is inverted by CELL2COL.
%   __________________________________________________________________
%
%   'col2cell --t' runs a test.
%
%   See also CELL2COL, COL2MAT, MAT2COL, COLBREAKS.
%
%   Usage: x=col2cell(col);
%          [x1,x2,...,xN]=col2cell(c1,c2,...,cN);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2014 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    col2cell_test,return
end

i=0;
if ~isempty(varargin{1})
    while i<nargin
        i=i+1;
        col=varargin{i};
        if i==1
            %Note minor changes to make this work on matrices rather than
            %just columns.  Not currently inverted by CELL2COL.
            ia=[1;find(isnan(col(1:end-1,1)))+1];
            ib=find(isnan(col(:,1)))-1;
            %[L,ia,ib]=blocklen(col);
            ray=cell(length(ia),1);
        end
        for j=1:length(ia)
            if ~isempty(col)
                ray{j}=col(ia(j):ib(j),:);
            else
                ray{j}=[];
            end
        end
        varargout{i}=ray;
    end
    if nargout==0 
        eval(to_overwrite(nargin));
     end
else
    for i=1:nargout
        varargout{i}=varargin{1};
    end
end


 

function[]=col2cell_test
c1=   [0 0 0 nan 1 1 nan 3 nan 2 2 2 nan]';
c2=   [3 3 3 nan 2 2 nan nan nan 3 nan 3 nan]';
c3=   [1 1 1 nan 2 2 nan 3 nan 4 4 4 nan]';
[d1,d2,d3]=col2cell(c1,c2,c3);
[e1,e2,e3]=cell2col(d1,d2,d3);
reporttest('COL2CELL inverted by CELL2COL',aresame(c1,e1)&&aresame(c3,e3))
%reporttest('COL2CELL inverted by CELL2COL',aresame(c1,e1)&&aresame(c2,e2)&&aresame(c3,e3))

col2cell(c1,c2,c3);
cell2col(c1,c2,c3);
reporttest('COL2CELL inverted by CELL2COL with overwriting',aresame(c1,e1)&&aresame(c2,e2)&&aresame(c3,e3))

c1=   [0 0 0 nan 1 1 nan 3 nan 2 2 2 nan]';
c2=   [3 3 3 nan 2 2 nan nan nan 3 nan 3 nan]';
c3=   [1 1 1 nan 2 2 nan 3 nan 4 4 4 nan]';
c1=c1+sqrt(-1)*c1/2;
c2=c2+sqrt(-1)*c2*2;
c3=c3+sqrt(-1)*c3/3;

[d1,d2,d3]=col2cell(c1,c2,c3);
[e1,e2,e3]=cell2col(d1,d2,d3);
reporttest('COL2CELL inverted by CELL2COL for complex data',aresame(c1,e1)&&aresame(c3,e3))
%reporttest('COL2CELL inverted by CELL2COL',aresame(c1,e1)&&aresame(c2,e2)&&aresame(c3,e3))

col2cell(c1,c2,c3);
cell2col(c1,c2,c3);
reporttest('COL2CELL inverted by CELL2COL for complex data with overwriting',aresame(c1,e1)&&aresame(c2,e2)&&aresame(c3,e3))

d1=col2cell([c1 c1]);

