list=dir('C:\Users\stationcalcul\GoogleDrive\MEOP-CTD\final_dataset_hr\prof\*.nc');
I_all=[];
J_all=[];
temp_HR_LR=[];
sal_HR_LR=[];
dens_HR_LR=[];
diff_temp_LR=[];
diff_sal_LR=[];
diff_dens_LR=[];
diff_temp_HR=[];
diff_sal_HR=[];
diff_dens_HR=[];

for ll=1:length(list)
    nc_nameHR = ['C:\Users\stationcalcul\GoogleDrive\MEOP-CTD\final_dataset_hr\prof\' list(ll).name];
    nc_nameLR=char(regexpdir('C:\Users\stationcalcul\GoogleDrive\MEOP-CTD\final_dataset_ncARGO\',list(ll).name(1:end-14)));
    Procedure
end


subplot(2,3,1)
plot(nanmean(temp_HR_LR'),-p,'-k',nanmean(temp_HR_LR(:,I_all)'),-p,'-b',nanmean(temp_HR_LR(:,J_all)'),-p,'-r')
xlabel('Mean diff temp')
ylabel('Depth')
subplot(2,3,2)
plot(nanmean(sal_HR_LR'),-p,'-k',nanmean(sal_HR_LR(:,I_all)'),-p,'-b',nanmean(sal_HR_LR(:,J_all)'),-p,'-r')
xlabel('Mean diff sal')
subplot(2,3,3)
plot(nanmean(dens_HR_LR'),-p,'-k',nanmean(dens_HR_LR(:,I_all)'),-p,'-b',nanmean(dens_HR_LR(:,J_all)'),-p,'-r')
legend('all','LR -10pts','LR +10pts','Location','southwest','FontSize',2)
xlabel('Mean diff sigma0')
subplot(2,3,4)
plot(nanstd(temp_HR_LR'),-p,'-k',nanstd(temp_HR_LR(:,I_all)'),-p,'-b',nanstd(temp_HR_LR(:,J_all)'),-p,'-r')
xlabel('STD diff temp')
ylabel('Depth')
subplot(2,3,5)
plot(nanstd(sal_HR_LR'),-p,'-k',nanstd(sal_HR_LR(:,I_all)'),-p,'-b',nanstd(sal_HR_LR(:,J_all)'),-p,'-r')
xlabel('STD diff temp')
subplot(2,3,6)
plot(nanstd(dens_HR_LR'),-p,'-k',nanstd(dens_HR_LR(:,I_all)'),-p,'-b',nanstd(dens_HR_LR(:,J_all)'),-p,'-r')
xlabel('STD diff temp')
suptitle(['individu : ' nc_att.smru_platform_code ' Profil mean and std difference HR-LR with number points LR threshold = 10']);
print([dirplot '/mean_std_diff_nbpts_LR_all'],'-depsc');
print([dirplot '/mean_std_diff_nbpts_LR_all'],'-dpng','-r400');


%%
subplot(2,3,1)
plot(nanmean(diff_temp_HR'),-p,'-g',nanmean(diff_temp_LR'),-p,'-k',nanmean(diff_temp_LR(:,I_all)'),-p,'-b',nanmean(diff_temp_LR(:,J_all)'),-p,'-r')
xlabel('Mean diff temp')
ylabel('Depth')
subplot(2,3,2)
plot(nanmean(diff_sal_HR'),-p,'-g',nanmean(diff_sal_LR'),-p,'-k',nanmean(diff_sal_LR(:,I_all)'),-p,'-b',nanmean(diff_sal_LR(:,J_all)'),-p,'-r')
xlabel('Mean diff sal')
subplot(2,3,3)
plot(nanmean(diff_dens_HR'),-p,'-g',nanmean(diff_dens_LR'),-p,'-k',nanmean(diff_dens_LR(:,I_all)'),-p,'-b',nanmean(diff_dens_LR(:,J_all)'),-p,'-r')
legend('HR','all LR','LR -10pts','LR +10pts','Location','southwest','FontSize',2)
xlabel('Mean diff sigma0')
subplot(2,3,4)
plot(nanstd(diff_temp_HR'),-p,'-g',nanstd(diff_temp_LR'),-p,'-k',nanstd(diff_temp_LR(:,I_all)'),-p,'-b',nanstd(diff_temp_LR(:,J_all)'),-p,'-r')
xlabel('STD diff temp')
ylabel('Depth')
subplot(2,3,5)
plot(nanstd(diff_sal_HR'),-p,'-g',nanstd(diff_sal_LR'),-p,'-k',nanstd(diff_sal_LR(:,I_all)'),-p,'-b',nanstd(diff_sal_LR(:,J_all)'),-p,'-r')
xlabel('STD diff temp')
subplot(2,3,6)
plot(nanstd(diff_dens_HR'),-p,'-g',nanstd(diff_dens_LR'),-p,'-k',nanstd(diff_dens_LR(:,I_all)'),-p,'-b',nanstd(diff_dens_LR(:,J_all)'),-p,'-r')
xlabel('STD diff temp')
suptitle(['all individu : Profil mean and std difference cor - raw']);
print([dirplot '/mean_std_diff_cor-raw_all'],'-depsc');
print([dirplot '/mean_std_diff_cor-raw_all'],'-dpng','-r400');
