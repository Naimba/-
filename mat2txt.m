function [] = mat2txt(ch,X)
%ch:保存文件的文件路径
%X:保存的数组名称

fid=fopen(ch,'wt'); %写入的文件，各函数后面有说明
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
disp(['保存完毕,路径在',ch])