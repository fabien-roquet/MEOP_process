function generate_plot2(conf,EXP,one_smru_name, suffix, mode)

if isempty(conf),
    conf = init_mirounga;
end

if ~exist('one_smru_name','var') % all tags from EXP deployment
    one_smru_name = '';
elseif isempty(EXP),
    EXP=EXP_from_smru_name(one_smru_name);
end

close all

if ~exist('one_smru_name','var') % all tags from EXP deployment
    info_deployment=load_info_deployment(conf,EXP);
    if isempty(info_deployment.list_smru_name)
        return
    end
else  % tag smru_tag only
    info_deployment=load_info_deployment(conf,EXP,one_smru_name);
    if isempty(info_deployment.list_smru_name)
        return
    end
end


%% load configuration of plots
do_figure='on';
plot_conf=[];
if ~any(strcmp(conf.table_param.Properties.RowNames,EXP)),
    plot_conf.temp_error=0.1;
    plot_conf.psal_error=0.2;
else
    for field = conf.table_param.Properties.VariableNames
        plot_conf = setfield(plot_conf,field{1},conf.table_param{EXP,field{1}});
    end
end

if isfield(plot_conf,'Tcontour') 
    if isempty(plot_conf.Tcontour{1})
        plot_conf = rmfield(plot_conf,'Tcontour');
    else
        plot_conf.Tcontour = eval(plot_conf.Tcontour{1});
    end
end

if isfield(plot_conf,'Scontour') 
    if isempty(plot_conf.Scontour{1})
        plot_conf = rmfield(plot_conf,'Scontour');
    else
        plot_conf.Scontour = eval(plot_conf.Scontour{1});
    end
end

fields = fieldnames(plot_conf)';
for field = fields
    if any(isnan(getfield(plot_conf,field{1})))
        plot_conf = rmfield(plot_conf,field{1});
    end
end


%% do the plots
if exist('suffix','var') & exist('mode','var'),

    disp(['plot diag ' EXP ': ' suffix(2:end) ', ' mode]);
    sc_plot_data_tags;
    sc_plot_data_deploy;
    sc_build_latex;

else

    if length(info_deployment.list_tag_lr0)
        suffix = '_lr0'; mode = 'raw';
        disp(['plot diag ' EXP ': ' suffix(2:end) ', ' mode]);
        sc_plot_data_tags;
        sc_plot_data_deploy;
        sc_build_latex;
    end

    if length(info_deployment.list_tag_lr1)
        suffix = '_lr1'; mode = 'adj';
        disp(['plot diag ' EXP ': ' suffix(2:end) ', ' mode]);
        sc_plot_data_tags;
        sc_plot_data_deploy;
        sc_build_latex;
    end

    if length(info_deployment.list_tag_lr1)
        suffix = '_hr1'; mode = 'adj';
        disp(['plot diag ' EXP ': ' suffix(2:end) ', ' mode]);
        sc_plot_data_tags;
        sc_plot_data_deploy;
        sc_build_latex;
    end

    if length(info_deployment.list_tag_hr2)
        suffix = '_hr2'; mode = 'adj';
        disp(['plot diag ' EXP ': ' suffix(2:end) ', ' mode]);
        sc_plot_data_tags;
        sc_plot_data_deploy;
        sc_build_latex;
    end

    if length(info_deployment.list_tag_fr1)
        suffix = '_fr1'; mode = 'adj';
        disp(['plot diag ' EXP ': ' suffix(2:end) ', ' mode]);
        sc_plot_data_tags;
        sc_plot_data_deploy;
        sc_build_latex;
    end

end

