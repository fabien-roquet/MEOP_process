% PRES, TEMP, SALI, CHLA, DOXY, LIGHT, std_std --> **_std_lev

% std_1dbar = (0:2000)';
N_std = length(std_lev);
N_prof = size(PRES,2);

PRES_std_lev = repmat(std_lev,1,N_prof);
TEMP_std_lev = NaN*PRES_std_lev;
SALI_std_lev = NaN*PRES_std_lev;
CHLA_std_lev = NaN*PRES_std_lev;
DOXY_std_lev = NaN*PRES_std_lev;
LIGHT_std_lev = NaN*PRES_std_lev;

% [Sgrid,Tgrid] = meshgrid(10:50,-4:.5:25);
% Dgrid = sw_dens0(Sgrid,Tgrid);

for kprof=1:N_prof,
    
    z=PRES(:,kprof);
    t=TEMP(:,kprof);
    s=SALI(:,kprof);
    f=CHLA(:,kprof);
    o=DOXY(:,kprof);
    l=LIGHT(:,kprof);
        
    d = sw_dens0(s,t);
    t_std = NaN*zeros(N_std,1);
    s_std = NaN*zeros(N_std,1);
    f_std = NaN*zeros(N_std,1);
    o_std = NaN*zeros(N_std,1);
    l_std = NaN*zeros(N_std,1);
%     t_1dbar=[];
%     s_1dbar=[];

    t_in_T=[];
    if ~isempty(t)&length(find(~isnan(z.*t)))>1
        z_in_T=z(~isnan(t) & ~isnan(z));
        t_in_T=t(~isnan(t) & ~isnan(z));
        [z_in_T,I]=sort(z_in_T);
        t_in_T=t_in_T(I);
        I=[find(z_in_T(1:end-1)~=z_in_T(2:end));length(z_in_T)];
        z_in_T=z_in_T(I);
        t_in_T=t_in_T(I);
        if length(t_in_T)>5 %...expected to avoid isolated values
            t_std = interp1(z_in_T,t_in_T,std_lev);
            % t_1dbar = interp1(z_in,t_in,std_1dbar);
        end
    end
    
    if ~isempty(s)&length(find(~isnan(z.*s)))>1
        z_in=z(~isnan(s) & ~isnan(z));
        s_in=s(~isnan(s) & ~isnan(z));
        [z_in,I]=sort(z_in);
        s_in=s_in(I);
        I=[find(z_in(1:end-1)~=z_in(2:end));length(z_in)];
        z_in=z_in(I);
        s_in=s_in(I);
        if length(s_in)>5 %...expected to avoid isolated values
            s_std = interp1(z_in,s_in,std_lev);
            % use gsw_SA_CT_interp for S and T instead of linear interp if possible
            if length(t_in_T)>5
                z_in_TS = sort(unique([z_in_T;z_in]));
                t_in_TS = interp1(std_lev,t_std,z_in_TS);
                s_in_TS = interp1(std_lev,s_std,z_in_TS);
                [s_std, t_std] = gsw_SA_CT_interp(s_in_TS,t_in_TS,z_in_TS,std_lev);
                % s_1dbar = interp1(z_in,s_in,std_1dbar);
            end
            
        end
    end
        
    if ~isempty(f)&length(find(~isnan(z.*f)))>1
        z_in=z(~isnan(f) & ~isnan(z));
        f_in=f(~isnan(f) & ~isnan(z));
        [z_in,I]=sort(z_in);
        f_in=f_in(I);
        I=[find(z_in(1:end-1)~=z_in(2:end));length(z_in)];
        z_in=z_in(I);
        f_in=f_in(I);
        if length(f_in)>5 %...expected to avoid isolated values
            f_std = interp1(z_in,f_in,std_lev);
            % f_1dbar = interp1(z_in,f_in,std_1dbar);
        end
    end
    
    if ~isempty(o)&length(find(~isnan(z.*o)))>1
        z_in=z(~isnan(o) & ~isnan(z));
        o_in=o(~isnan(o) & ~isnan(z));
        [z_in,I]=sort(z_in);
        o_in=o_in(I);
        I=[find(z_in(1:end-1)~=z_in(2:end));length(z_in)];
        z_in=z_in(I);
        o_in=o_in(I);
        if length(o_in)>5 %...expected to avoid isolated values
            o_std = interp1(z_in,o_in,std_lev);
            % o_1dbar = interp1(z_in,o_in,std_1dbar);
        end
    end
    
    if ~isempty(l)&length(find(~isnan(z.*l)))>1
        z_in=z(~isnan(l) & ~isnan(z));
        l_in=l(~isnan(l) & ~isnan(z));
        [z_in,I]=sort(z_in);
        l_in=l_in(I);
        I=[find(z_in(1:end-1)~=z_in(2:end));length(z_in)];
        z_in=z_in(I);
        l_in=l_in(I);
        if length(l_in)>5 %...expected to avoid isolated values
            l_std = interp1(z_in,l_in,std_lev);
            % o_1dbar = interp1(z_in,o_in,std_1dbar);
        end
    end
    
    TEMP_std_lev(:,kprof) = t_std;
    SALI_std_lev(:,kprof) = s_std;
    CHLA_std_lev(:,kprof) = f_std;
    DOXY_std_lev(:,kprof) = o_std;
    LIGHT_std_lev(:,kprof) = l_std;
end

