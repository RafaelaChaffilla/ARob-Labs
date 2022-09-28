%Guarda os dados das variáveis do simulink em ficheiros de texto

t = height.time;
h = height.signals.values;
p = pic.signals.values; 
fileID=fopen('dados3.txt','w');
for i = 1:length(t)
    fprintf(fileID,'%f   %f   %f\n', t(i),h(i), p(i));
end
fclose(fileID);