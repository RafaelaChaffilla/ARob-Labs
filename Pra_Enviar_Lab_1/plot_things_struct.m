clear all;
close all;

%NOTA: para este c�digo funcionar:
%
%      Nome dos ficheiros: 
%
%      1 -  devem come�ar com a indica��o do k e do h, e deve estar
%      exatamente da mesma forma escrito para o experimental e o teorico
%      
%      2 -  depois do k e h, dizer "_teorico" ou "_experimental"
%      OBRIGAT�RIO

path(pathdef);
path = [pwd '\Dados\Theta\'];

i=1;
previous = 'lol';

files = dir(strcat(path,'*.mat'));
for file = files'
    dados(i).nome   = file.name;
    dados(i).KH     = erase(file.name,"_commands.mat");
    dados(i).KH     = erase(dados(i).KH,".mat");
    dados(i).KH     = strrep(dados(i).KH,'_',' ');
    holder          = load(strcat(path,file.name));
    dados(i).tempo  = holder.ans(1,:);
    if(size(holder.ans,1) == 2)
        dados(i).theta  = holder.ans(2,:);
    else
        dados(i).theta  = holder.ans(3,:);
    end
    i = i+1;
end

for coisa = dados
    if(strcmp(previous, coisa.KH) == 1) % o anterior tbm era os mesmo dados
        plot(coisa.tempo,coisa.theta);
        legend('experimental', 'comando', 'Location', 'southeast');
    else
        figure();
        plot(coisa.tempo,coisa.theta);
        hold on;
        title(coisa.KH);
        %hold off;
    end
    previous = coisa.KH;
    
end