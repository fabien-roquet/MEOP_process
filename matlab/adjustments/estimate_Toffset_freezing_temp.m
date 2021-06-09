function Toffset = estimate_Toffset_freezing_temp(argo_qc,nosal)



argo_qc.FREEZING_LINE = sw_fp(argo_qc.PSAL,0);
if nosal
    argo_qc.FREEZING_LINE = argo_qc.TEMP * 0 + sw_fp(34.5,0);
end
argo_qc.PTEMP = sw_ptmp(argo_qc.PSAL,argo_qc.TEMP,argo_qc.PRES,0);
aux = sw_ptmp(argo_qc.TEMP*0+34.5,argo_qc.TEMP,argo_qc.PRES,0);
I=find(isnan(argo_qc.PSAL)); argo_qc.PTEMP(I) = aux(I);
Toffset=zeros(argo_qc.ntag,1);

if length(argo_qc.JULD)>20,
    [m,i]=sort(argo_qc.PTEMP(2,:)-argo_qc.FREEZING_LINE(2,:));
    Toffset=ceil(median(m(1:10))*100)/100;
    if isnan(Toffset), Toffset=0; end
end

disp(sprintf('%5.2f', Toffset));





