
function [a]=readAbed(csvfile)   
%
% ��������� a  �����ݷֽ�� 5 �У��������Ե��Լ����ܶ� a(:,1:5)
% ���룺csvfile  abed ����������ļ� 
%
% �������� 1
%   f='C:\Users\Think\Desktop\abed-1pop5replicator.csv'
%   [a data]=readAbed(csvfile)
%   a ���������ʱ������
%   [num,txt,raw] = xlsread('abed-1pop4d.csv');
%
%

    [num,txt,raw] = xlsread(csvfile);
    v = cell2table(raw(length(txt(:,1))+1:length(raw(:,1)),:));
    tm=table2array(v);
    data=v;
% �����ݷֽ�� 5 �У��������Ե��Լ����ܶ� a(:,1:5)
    N=4; 
        a=[];   
        for i =2:4:2+N*4-1
                a=[a tm(:,i)-tm(:,i+4)]; 
        end
        a=[a tm(:,2+N*4)];
% ���� 6 7 ����Ϊ�ϳɵĲ��ԣ�        
%         a(:,6) =  a(:,1) +  a(:,4);
%         a(:,7) =  a(:,2) +  a(:,5); 
end
