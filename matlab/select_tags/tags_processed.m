function EXPs = tags_processed(varargin);
%% choose list tag to be processed
% [Itag]=tags_processed(varargin);
% Arguments:
%   1) conf: required struct variable obtained when initializing sc_mirounga
%   2) if cell: list of deployments that need to be processed
%      if str: mode of selection ('unprocessed_tags' or 'new_tags' or a valid NATION)
% Output:
%   1) EXPs: list of deployments that must be processed
%
% tag_processed(): list info for all processed tags
% tag_processed(one_smru_name): info on one_smru_name
% tag_processed(list_smru_name): info on tags in list_smru_name (cell string array)
% tag_processed(conf,list_smru_name): info on tags in list_smru_name
% tag_processed(conf,'unprocessed_tags'): list info on unprocessed_tags
% tag_processed(conf,'new_tags'): list info on unprocessed_tags
% tag_processed(conf,'NATION',NATION): list info on tags from NATION
%
% conf is returned by function init_mirounga. Autoimatically retrieved if conf=[].

Itag=[];
if nargin==0
    conf = init_mirounga;
    EXPs = conf.list_deployment(conf.list_deployment.process==1,:);
    EXPs = sortrows(EXPs,{'country','pi_code','start_date_jul'});
    
elseif nargin==1
    smru_name = varargin{1};
    conf = init_mirounga;
    if isstr(smru_name),
        smru_name = {smru_name};
    end
    EXPs = conf.list_deployment(smru_name,:);
    
elseif nargin==2 && iscell(varargin{2})
    conf = varargin{1};
    if isempty(conf), conf=init_mirounga; end
    list = varargin{2};
    EXPs = conf.list_deployment(list,:);
    
elseif nargin==2 && isstr(varargin{2})
    conf = varargin{1};
    if isempty(conf), conf=init_mirounga; end
    mode = varargin{2};
    switch mode
        case 'unprocessed_tags' % tags that have not been processed
            
            EXPs = tags_processed(conf);
            count = 0;
            list = {};
            for kEXP = 1:length(EXPs.deployment_code),
                EXP = EXPs.deployment_code{kEXP};
                info_deployment = load_info_deployment(conf,EXP);
                if exist([conf.datadir,info_deployment.NATION,'/',info_deployment.EXP],'file'),
                    count = count+1;
                    list{count} = EXP;
                end
            end
            EXPs(list,:)=[];
            
        case 'new_tags' % tags that were not present in previous release
            
            EXPs = conf.list_deployment(conf.list_deployment.process==1 & ...
                strcmp(conf.list_deployment.last_version,''),:);

        otherwise % invalid mode
            
            error('Invalid selection mode in tag_processed');
            
    end
    EXPs=sortrows(EXPs,{'country','pi_code','start_date_jul'});

elseif nargin==3 && isstr(varargin{2}) && strcmp(varargin{2},'NATION')
    conf = varargin{1};
    if isempty(conf), conf=init_mirounga; end
    NATION = varargin{3};
    EXPs = conf.list_deployment(conf.list_deployment.process==1 & ...
        strcmp(conf.list_deployment.country,NATION),:);
    
else
    error('Too many arguments in tag_processed');
    
end



