syms a b;
rM5ab=[ a 1 2 3 4+b
        4 a 1 2+b 3
        3 4 a+b 1 2
        2 3+b 4 a 1
        1+b 2 3 4 a];  
eva = [];eve = []; Mabed=[];
gen_abed_payoff = []; 
for k=2:0.01:3; 
        a = k; b=0; 
        eval_M5ab = eval(rM5ab)
%         [eigen_vector, d_eigen_value, rest_point_no] = ...
%             get_eigensystem_from_5x5(eval_M5ab); 
%         eva = [eva d_eigen_value];
%         eve = [eve; eigen_vector]; 
%               str1 = [98 98 eval_M5ab(1,:) 99 98 eval_M5ab(2,:)  99 98 eval_M5ab(3,:)  99 98 eval_M5ab(4,:)  99 98 eval_M5ab(5,:) 99 99];
%         Mabed = [Mabed; str1];
              str2 = strcat(' "[[', num2str(eval_M5ab(1,:)), '] [', num2str(eval_M5ab(2,:)), '] [', num2str(eval_M5ab(3,:)), '] [',num2str(eval_M5ab(4,:)), '] [',num2str(eval_M5ab(5,:)), ']]" ');
        gen_abed_payoff = [gen_abed_payoff str2];  
end
gen_abed_payoffu =  strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(gen_abed_payoff,'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' ') ;
gen_abed_payoffu

% netlogo abed InvarMani_diag100 

syms a b;
rM5ab=[ a 1+b 2 3 4
        4 a 1+b 2 3
        3 4 a 1+b 2
        2 3 4 a 1+b
        1+b 2 3 4 a];  
eva = [];eve = []; Mabed=[];
gen_abed_payoff = []; 
for k=2:0.5:3; for b=-1:1:1; 
        a = k; 
        eval_M5ab = eval(rM5ab)
        [eigen_vector, d_eigen_value, rest_point_no] = ...
            get_eigensystem_from_5x5(eval_M5ab); 
        eva = [eva d_eigen_value];
        eve = [eve; eigen_vector]; 
              str1 = [98 98 eval_M5ab(1,:) 99 98 eval_M5ab(2,:)  99 98 eval_M5ab(3,:)  99 98 eval_M5ab(4,:)  99 98 eval_M5ab(5,:) 99 99];
        Mabed = [Mabed; str1];
              str2 = strcat(' "[[', num2str(eval_M5ab(1,:)), '] [', num2str(eval_M5ab(2,:)), '] [', num2str(eval_M5ab(3,:)), '] [',num2str(eval_M5ab(4,:)), '] [',num2str(eval_M5ab(5,:)), ']]" ');
        gen_abed_payoff = [gen_abed_payoff str2];  
    end;end
gen_abed_payoffu =  strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(gen_abed_payoff,'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' ') ;
gen_abed_payoffu


% ab121
 syms a b
 gen_abed_payoff=[]; gameid=1;
rM5ab=[ a-b 1+b 2 3 4
        4 a-b 1+b 2 3
        3 4 a-b 1+b 2
        2 3 4 a-b 1+b
        1+b 2 3 4 a-b];
    
    for a=0:0.5:5; for b=0:0.5:5; 
            
            eval_M5ab = eval(rM5ab)  
                 str2 = strcat(' "[[', num2str(eval_M5ab(1,:)), '] [', num2str(eval_M5ab(2,:)), '] [', num2str(eval_M5ab(3,:)), '] [',num2str(eval_M5ab(4,:)), '] [',num2str(eval_M5ab(5,:)), ']]" ');
           gen_abed_payoff = [gen_abed_payoff str2];
           Mab121(:,:,gameid) = eval_M5ab; 
           gameid=gameid+1;
    end;end
save('Mab121.mat','Mab121')
gen_abed_payoffu =  strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(gen_abed_payoff,'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' ') ;
gen_abed_payoffu
 

% double_center_manifold 
 syms a b
 gen_abed_payoff=[]; gameid=1;
 rM5ab=[0 a 1 -1 -a;
        -a 0 a 1 -1;
        -1 -a 0 a 1;
        1 -1 -a 0 a;
        a 1 -1 -a 0];
    for a=-5:0.1:5            
            eval_M5ab = eval(rM5ab)  
                 str2 = strcat(' "[[', num2str(eval_M5ab(1,:)), '] [', num2str(eval_M5ab(2,:)), '] [', num2str(eval_M5ab(3,:)), '] [',num2str(eval_M5ab(4,:)), '] [',num2str(eval_M5ab(5,:)), ']]" ');
           gen_abed_payoff = [gen_abed_payoff str2];
           Mab121(:,:,gameid) = eval_M5ab; 
           gameid=gameid+1;
    end; 
double_center_manifold_M121str =  strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(gen_abed_payoff,'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' ') ;
double_center_manifold_M121str
double_center_manifold_M121matrix = Mab121; 
save('F:\tmpdata\parameter202107\double_center_manifold_payoff.mat','double_center_manifold_M121str','double_center_manifold_M121matrix')

function fj_in_order =  fj_in_order(Tm)
load('F:\tmpdata\parameter202107\double_center_manifold_payoff.mat')

index=0;
for j=1:size(double_center_manifold_M121matrix,3)
    if sum(sum(abs(double_center_manifold_M121matrix(:,:,j) - Tm)))<0.01; index=j;end
end 
fj_in_order = index;
end








