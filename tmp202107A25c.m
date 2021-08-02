 load('F:\tmpdata\result_ana202107\InvMani0725t001_200_onlyMean.mat')
a=[];for k=1:200; a=[a;rr_E_mean_std(:,1,k)'];end
a(:,6:7) =[min(a')' max(a')'];
xx=find(a(:,6) > 0.10)
      for kk=1:size(xx)
          fj=xx(kk)
          E_angularmom = rr_E_angularmom(:,:,fj);
      % theoretical result
         [eigen_vector, d_eigen_value, T_eigencycle, n0, solu_list] = ...
                 get_eigensystem_from_5x5(rr_M5x5(:,:,fj)); 
       % statics analysis
         [coef_coefint3col,resid_residint3col,R2_F_p_evar] = ... 
                 multi_reg_Eam_Tec(E_angularmom,T_eigencycle);   
                  bbb=coef_coefint3col'  
                  R2_F_p_evar
      end