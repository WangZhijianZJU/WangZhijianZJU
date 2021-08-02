
load('F:\tmpdata\result_ana202107\InvMani0725_onlyMean.mat')
mean_std10 = []; for k=1:489; mean_std10 = [mean_std10; rr_E_mean_std(:,1,k)']; end
mean_std10(:,6:7) =[min(mean_std10')' max(mean_std10')'];
% 
bbb=[]; singularity =[]; eigen_value_ri=[];
xx=find(mean_std10(:,6) > 0.12)
for kk=1:size(xx)
    fj=xx(kk)
    E_angularmom = rr_E_angularmom(:,:,fj);
    % theoretical result
    [eigen_vector, d_eigen_value, T_eigencycle, n0, solu_list] = ...
    get_eigensystem_from_5x5(rr_M5x5(:,:,fj));
    % statics analysis 
    [coef_coefint3col,resid_residint3col,R2_F_p_evar] = ...
    multi_reg_Eam_Tec(E_angularmom,T_eigencycle);
    if size(coef_coefint3col,1) == 6; 
        eigen_value_ri=[real(d_eigen_value)' imag(d_eigen_value)'];
        bbb=[bbb; mean_std10(fj,:) 99 coef_coefint3col(:,1)' 99 R2_F_p_evar 97 eigen_value_ri 96 norm(E_angularmom) fj ]
    else
        singularity =[singularity; fj];
    end
end
save('F:\tmpdata\result_ana202107\InvMani0725_onlyMeanbbb12.mat','bbb','singularity' ) 


% xxx = bbb(:,22:25)
% xx1 = xxx+abs(xxx);
% xx2=xx1/2; xx3=sqrt(xx2(:,1).*xx2(:,2))+ sqrt(xx2(:,3).*xx2(:,4))
% scatter(xx3, bbb(:,32))


% in=[];for j=1:size(xx3); tmp=find(abs(xxx(j,:)-xx3(j))<0.00001); if size(tmp)>0; tmp1=tmp(1);else;tmp1=0;end;in=[in; tmp1]; end
% yyy = bbb(:,27:30);
% yy1=zeros(size(yyy,1),1);
% for j=1:size(yyy,1); if in(j,1)>0; yy1(j)=yyy(j,in(j,1));end;end

y = rrr(:,1);X = [ones(size(y)) ,...
rrr(:,2),...
rrr(:,3) ];
[b,bint,r,rint,stats] = regress(y,X)