function f = latex2MxWithMxPrecision(A, precision)
% û���뾫�Ȳ���ʱ��Ĭ�Ͼ���ΪС����4λ
if nargin == 1
    precision = '4';
else
    precision = int2str(precision);
end
% ���嵥һԪ�������ʽ
out_num = [' %0.' precision 'f &'];
% ������������ж�
z = zeros(1, str2num(precision) + 1);
z(1) = '.';
z(2 : end) = '0';
z = char(z);
% ������С
[r c] = size(A);
nc = zeros(1, c);
nc(:) = 99;  % ���character c
% ���ɵ�һ��Latex���
out = sprintf('\\left(\n\t\\begin{array}{%s}', char(nc));
% ����ѭ���������������������Latex���
for i = 1 : r
    out = [out sprintf('\n\t')]; % ����
    for j = 1 : c
        temp = sprintf(out_num, A(i, j));
        % С��λ��Ϊ��ʱ������ȡ������1.0001ȡΪ1
        dot_position = find(temp == '.');
        if temp(dot_position : end - 2) == z
            temp = temp(1 : dot_position - 1);
            temp = [temp ' &'];
            % Ҫȡ��ʱ�����и��ţ�����붪��
%             if temp(2) == '-'
%                temp = [temp(1) temp(3 : end)]; 
%             end
        end
        out = [out temp];
    end
    % ��������'&'��
    out = out(1 : end - 1);
    % ��ĩ����'\\'��
    out = [out '\\'];
end
% �������һ���������
out = [out sprintf('\n\t\\end{array}\n\\right)')];
f = out;