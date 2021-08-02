
% Programmer: socexp, WZJ, 20210724 simplified, rewrite 
% Notice: Only for 5 strategy case
% Input: abedfilesname, full path and file name
% Output: density of state x at time t, 5 column  

function dos_x_t = get_dos_x_t_from_abeddata(abedfilesname)
            [num,txt,raw] = xlsread(abedfilesname);
                            t=length(txt(:,1));
                            u=raw(t+1:length(raw(:,1)),:);
                            v=cell2table(u);
                            ps1 = [v.u1 v.u2 v.u6 v.u10  v.u14 v.u18];
             dos_x_t = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18];
end

