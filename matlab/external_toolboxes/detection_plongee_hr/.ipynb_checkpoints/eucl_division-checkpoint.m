function[q,r]=eucl_division(a,b)
%division euclidienne de a par b
% a=b*q +r
%--------------------------------------------------------------------------
q=floor(a/b);
r=a-b*q;