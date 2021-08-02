
function [a]=readAbed(csvfile)   
%
% 输出：数组 a  将数据分解成 5 列，各个策略的自己的密度 a(:,1:5)
% 输入：csvfile  abed 的输出数据文件 
%
% 测试用例 1
%   f='C:\Users\Think\Desktop\abed-1pop5replicator.csv'
%   [a data]=readAbed(csvfile)
%   a 就是所需的时间序列
%   [num,txt,raw] = xlsread('abed-1pop4d.csv');
%
%

    [num,txt,raw] = xlsread(csvfile);
    v = cell2table(raw(length(txt(:,1))+1:length(raw(:,1)),:));
    tm=table2array(v);
    data=v;
% 将数据分解成 5 列，各个策略的自己的密度 a(:,1:5)
    N=4; 
        a=[];   
        for i =2:4:2+N*4-1
                a=[a tm(:,i)-tm(:,i+4)]; 
        end
        a=[a tm(:,2+N*4)];
% 构造 6 7 二列为合成的策略；        
%         a(:,6) =  a(:,1) +  a(:,4);
%         a(:,7) =  a(:,2) +  a(:,5); 
end
