function [] = mat2txt(ch,X)
%ch:�����ļ����ļ�·��
%X:�������������

fid=fopen(ch,'wt'); %д����ļ���������������˵��
[m,n]=size(X);
for i=1:1:m
    for j=1:1:n
       if j==n
         fprintf(fid,'%g\n',X(i,j));
      else
        fprintf(fid,'%g\t',X(i,j));
       end
    end
end
fclose(fid);
disp(['�������,·����',ch])