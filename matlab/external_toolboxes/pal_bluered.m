function pal=pal_bluered(N)

if ~exist('N','var')
    N=100;
end
%tab=[0 0 1 1;30 .2 .2 1; 49 1 1 1; 51 1 1 1; 70 1 .2 .2; 100 1 1 0];
tab=[0 0 1 1;30 .2 .2 1; 50 1 1 1; 70 1 .2 .2; 100 1 1 0];
pal=interp1(tab(:,1),tab(:,2:4),0:100/(N-1):100);

if mod(N,2)==1
    pal(ceil(N/2),:)=[1 1 1];
else
    pal(N/2,:)=[1 1 1];
    pal(N/2+1,:)=[1 1 1];
end
