% abed-v5-InvarMani.nlogo

load('manifold5try.mat'); 
gen_abed_payoff = [];  

%     for b=1:5:501 %size(R,1);  %202107242350  1-101
%     for b=506:5:1000 %size(R,1);  %202107250910  102-200
%     for b=1001:5:1500 %size(R,1);  %202107250910  201-300 
    for b=1501:5:2455 %size(R,1);  %202107250910  301-491   
        eval_M5ab = R(b:b+4,:); 
              str2 = strcat(' "[[', num2str(eval_M5ab(1,:)), '] ...
              [', num2str(eval_M5ab(2,:)), '] [', num2str(eval_M5ab(3,:)), '] ...
              [',num2str(eval_M5ab(4,:)), '] [',num2str(eval_M5ab(5,:)), ']]" ');
        gen_abed_payoff = [gen_abed_payoff str2]; 
    end;
 gen_abed_payoffu =  strrep(strrep(strrep(strrep(strrep(strrep(strrep( ...
     strrep(strrep(strrep(gen_abed_payoff,'  ',' '),'  ',' '),'  ',' '), ...
     '  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' '),'  ',' ') ;
gen_abed_payoffu


