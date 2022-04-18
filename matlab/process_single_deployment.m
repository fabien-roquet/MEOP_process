function process_single_deployment(EXP)
% process_single_deployment(EXP)
% run all the procedure for one tag only

conf = init_mirounga;
info_deployment=load_info_deployment(conf,EXP);
if isfield(info_deployment,'invalid_code') & info_deployment.invalid_code,
    return
end
import_ODV_data(EXP);    
remove_deployment(conf,EXP);
create_ncargo(conf,EXP);
create_fr0(conf,EXP);
create_fr0_without_lr0(conf,EXP);
update_metadata(conf,EXP);
apply_adjustments(conf,EXP);
apply_tlc(conf,EXP);
apply_tlc_fr(conf,EXP);
create_hr2(conf,EXP);
generate_odv4(conf,EXP);

generate_plot1(conf,EXP);
generate_plot2(conf,EXP);

visualize_tags(conf,EXP);