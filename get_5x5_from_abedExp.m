
%Programmer: socexp, ZJU, 20210723
% Input: experiment results from file_path
%Output: fet the payoff matrix from data as M5x5 
% Try with command 
% file_path='F:\tmpdata\result_exp202107';[M5x5, fileNames]=get_5x5_from_abedExp(file_path)
%

function [M5x5, fileNames]=get_5x5_from_abedExp(file_path)
% file_path='F:\tmpdata\result_exp202107';
fileFolder=fullfile(file_path);
dirOutput=dir(fullfile(fileFolder,'*.csv'));
fileNames={dirOutput.name}';

% M=zeros(5,5,:)
for fj = 1:size(fileNames,1);
        x = fileNames{fj}; 
        y= strrep(strrep(strrep(x,'[',''),']',' '),'.csv','');
        z=str2num(y);
        rMatrix=[]; 
        for k=0:4; 
            rMatrix=[rMatrix;z(1+k*5:5+k*5) ]; 
        end;
        M5x5(:,:,fj)=rMatrix;
end
