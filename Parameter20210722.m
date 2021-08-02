syms a b;
rM5ab=[ a 1 2 3 4+b
        4 a 1 2+b 3
        3 4 a+b 1 2
        2 3+b 4 a 1
        1+b 2 3 4 a];  
eva = [];eve = []; Mabed=[];
gen_abed_payoff = []; 
for k=2:0.5:3; 
    for b=-3:3:3; 
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
    end; 
end
gen_abed_payoffu =  strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(gen_abed_payoff,'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' ') ;
gen_abed_payoffu

