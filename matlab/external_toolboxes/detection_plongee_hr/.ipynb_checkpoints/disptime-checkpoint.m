function[time]=disptime(time)
%conversion d'une duree exprimee en secondes en heures, minutes et secondes
%--------------------------------------------------------------------------
[j,time]=eucl_division(time,3600*24);
[h,time]=eucl_division(time,3600);
[m,s]=eucl_division(time,60);
s=round(s);
time=[j,h,m,s];
