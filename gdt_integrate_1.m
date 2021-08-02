% Programmer: socexp, WZJ, 20210724 
% gdt: game dynamics theory
% Input: the subfold of the abed experiment csv files
% Output: statics analysis result 
%         save('manifold20210724try.mat', ...
%           'rr_abed_csv', ... % rr_abed_csv{:,:,3}
%         'rr_dos_x_t', ...
%         'rr_E_mean_std', ...
%         'rr_E_angularmom', ...
%         'rr_eigen_vector', ...
%         'rr_d_eigen_value', ...
%         'rr_T_eigencycle', ...
%         'rr_rest_no', ...
%         'rr_solu_list', .. %cell2mat(rr_solu_list(:,:,3))
%         'rr_coef_coefint3col', ...
%         'rr_resid_residint3col', ...
%         'rr_R2_F_p_evar', ...
%         'rr_M5x5');
%% Test command 

%gdt_integrate_1('F:\tmpdata\result_exp202107\','F:\tmpdata\result_ana202107\InvMani0725t001_200.mat')
% for ab121 case 20210737 
% gdt_integrate_1('F:\tmpdata\result_exp202107\ab121\','F:\tmpdata\result_ana202107\InvManiY65ab121.mat')
% gdt_integrate_1('F:\tmpdata\result_exp202107\double_center_manifold\','F:\tmpdata\result_ana202107\InvManiY65double_center_manifold.mat')
% for best 
% file_abedcsv_path = 'F:\tmpdata\result_exp202107\double_center_manifold\2centerbest-101\'; file_save_mat = 'F:\tmpdata\result_ana202107\InvManiY65double_center_best.mat' 
% for 2centerMu0.001-101
% file_abedcsv_path = 'F:\tmpdata\result_exp202107\double_center_manifold\2centerMu0.001-101\'; file_save_mat = 'F:\tmpdata\result_ana202107\2centerMu0.001-101.mat' 
%2centerMu0.10-101
% file_abedcsv_path = 'F:\tmpdata\result_exp202107\double_center_manifold\2centerMu0.10-101\'; file_save_mat = 'F:\tmpdata\result_ana202107\2centerMu0.10-101.mat' 

function csv101=gdt_integrate_1(file_abedcsv_path,file_save_mat)
            %file_abedcsv_path='F:\tmpdata\result_exp202107\';
         [rr_M5x5, fileNames]=get_5x5_from_abedExp(file_abedcsv_path); 
    
    for fj_in_abed = 1:size(fileNames,1);
        % abed experiment data 
            abed_csv = fileNames{fj_in_abed}; 
                abedfilesname = strcat(file_abedcsv_path,abed_csv)
            dos_x_t = get_dos_x_t_from_abeddata(abedfilesname);
         [E_angularmom,mn,am_eigencycleDim_t, E_mean, E_std] = ...
                from_N_colExp_out_am(dos_x_t);
                [E_mean; E_std] 
       % theoretical result
         [eigen_vector, d_eigen_value, T_eigencycle, n0, solu_list] = ...
                 get_eigensystem_from_5x5(rr_M5x5(:,:,fj_in_abed)); 
       % statics analysis
         [coef_coefint3col,resid_residint3col,R2_F_p_evar] = ... 
                 multi_reg_Eam_Tec(E_angularmom,T_eigencycle);   
                  coef_coefint3col'  
                  R2_F_p_evar
       % data organized
       fj =  fj_in_Torder(rr_M5x5(:,:,fj_in_abed));
          if size(d_eigen_value,1) == 5
            rr_abed_csv(:,:,fj) = {abed_csv}; % rr_abed_csv{:,:,3}
            rr_dos_x_t(:,:,fj) = dos_x_t;
            rr_E_mean_std(:,:,fj) = [E_mean' E_std'];
            rr_E_angularmom(:,:,fj) = E_angularmom;  
            rr_eigen_vector(:,:,fj) =  eigen_vector;
            rr_d_eigen_value(:,:,fj) = d_eigen_value'; 
            rr_T_eigencycle(:,:,fj) =  T_eigencycle;
            rr_rest_no(:,:,fj) =  n0;
            rr_solu_list(:,:,fj) = {solu_list}; %cell2mat(rr_solu_list(:,:,3))
            rr_coef_coefint3col(:,:,fj) = coef_coefint3col;
            rr_resid_residint3col(:,:,fj) = resid_residint3col;
            rr_R2_F_p_evar(:,:,fj) = R2_F_p_evar;
          end
    end 
    
    % save results
        %file_save_mat = 'manifold20210724try.mat'
        save(file_save_mat, 'rr_abed_csv', ...
                            'rr_dos_x_t', ...
                            'rr_E_mean_std', ...
                            'rr_E_angularmom', ...
                            'rr_eigen_vector', ...
                            'rr_d_eigen_value', ...
                            'rr_T_eigencycle', ...
                            'rr_rest_no', ...
                            'rr_solu_list', ...
                            'rr_coef_coefint3col', ...
                            'rr_resid_residint3col', ...
                            'rr_R2_F_p_evar', ...
                            'rr_M5x5'); 
    %report excel
    mean_l = []; for k=1:101; mean_l = [mean_l; rr_E_mean_std(:,1,k)']; end
    mean_l(:,6:7) =[min(mean_l')' max(mean_l')'];
    coef_eig = []; for k=1:101; coef_eig = [coef_eig; k rr_coef_coefint3col(:,1,k)' rr_R2_F_p_evar(1,:,k) 99 real(rr_d_eigen_value(1,:,k)) imag(rr_d_eigen_value(1,:,k)) norm(rr_E_angularmom(:,1,k))];end; 
    csv101=[mean_l coef_eig];
    g=csv101(:,[9:14 20:29]) ; 
    r2=[];for k=1:101; r1=[0];  for j=1:5; if abs(g(k,11+j))>0; r1=[r1 g(k,11+j)]; end; ;end; r2=[r2;r1];end; %try r1
    r2=[]; rr=[]; for k=1:101; r1=[0];  for j=1:5; if abs(g(k,11+j))>0; r1=[r1 g(k,11+j)]; else; rr=[rr; k j]; end; ;end; r2=[r2;r1];end; %try r1
    b3=[]; for k=1:101; r1=[g(k,1) 0];  for j=1:5; if abs(rr(k,2)-j)>0; r1=[r1 g(k,1+j)];end;end;b3=[b3;r1];end
    b3a=[b3(:,1:2)  b3(:,3)-b3(:,4) b3(:,1) b3(:,5)-b3(:,6)  b3(:,1)]; % b3(:,1) replace 0 ;
    b3a=[b3(:,1:2)  b3(:,3)-b3(:,4) b3(:,2) b3(:,5)-b3(:,6)  b3(:,2)]; % b3(:,2) replace 0 ;
    b_p_Im_Lnorm = [b3a csv101(:,15:24) r2 csv101(:,30) csv101(:,8)]; %csv(:,15:24) = id + 4p +5real + L;
    save(strcat(file_abedcsv_path,'b_p_Im_LnormMu0100.mat'),'b_p_Im_Lnorm')
end 
 

function rfj_in_Torder =  fj_in_Torder(M_in_abed)
load('F:\tmpdata\parameter202107\double_center_manifold_payoff.mat')

index=0;
for j=1:size(double_center_manifold_M121matrix,3)
    if sum(sum(abs(double_center_manifold_M121matrix(:,:,j) - M_in_abed)))<0.01; index=j;end
end 
rfj_in_Torder = index;
end
 
% load('F:\tmpdata\result_ana202107\InvManiY65double_center_manifold.mat'); 
% mean_std10 = []; for k=1:101; mean_std10 = [mean_std10; rr_E_mean_std(:,1,k)']; end
% mean_std10(:,6:7) =[min(mean_std10')' max(mean_std10')'];
% coef_6 = []; for k=1:101; coef_6 = [coef_6; k rr_coef_coefint3col(:,1,k)' rr_R2_F_p_evar(1,:,k) real(rr_d_eigen_value(1,:,k)) imag(rr_d_eigen_value(1,:,k)) norm(rr_E_angularmom(:,1,k))];end; 
% csv101=[mean_std10 coef_6];