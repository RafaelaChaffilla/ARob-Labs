clear all;
close all;

%NOTA: para este código funcionar:
%
%      Nome dos ficheiros: 
%
%      1 -  devem começar com a indicação do k e do h, e deve estar
%      exatamente da mesma forma escrito para o experimental e o teorico
%      
%      2 -  depois do k e h, dizer "_teorico" ou "_experimental"
%      OBRIGATÓRIO

path(pathdef);
path = [pwd '\Dados\'];

i=1;
previous = 'lol';

files = dir(strcat(path,'*.mat'));
for file = files'
    dados(i).nome = file.name;
    dados(i).KH = erase(file.name,"_teorico.mat");
    dados(i).KH = erase(dados(i).KH,"_experimental.mat");
    dados(i).KH = strrep(dados(i).KH,'_',' ');
    dados(i).dado = load(strcat(path,file.name));
    i = i+1;
end

for coisa = dados
    if(strcmp(previous, coisa.KH) == 1) % o anterior tbm era os mesmo dados
        plot(coisa.dado.ans(1,:),coisa.dado.ans(2,:));
        legend('experimental', 'teorico', 'Location', 'southeast');
    else
        figure();
        plot(coisa.dado.ans(1,:),coisa.dado.ans(2,:));
        hold on;
        title(coisa.KH);
        %hold off;
    end
    previous = coisa.KH;
    
end