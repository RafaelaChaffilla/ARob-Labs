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
path_1 = [pwd '\Dados\H\Raw\'];
path_2 = [pwd '\Dados\H\'];
addpath([pwd '\Functions']);

i=1;
previous = 'lol';

%% Carregar dados
files = dir(strcat(path_1,'*.mat'));
for file = files'
    dados(i).nome   = file.name;
    dados(i).KH     = erase(file.name,".mat");
    dados(i).KH     = erase(dados(i).KH,".mat");
    dados(i).KH     = strrep(dados(i).KH,'_',' ');
    holder          = load(strcat(path_1,file.name));

    dados(i).tempo  = holder.ans(1,:);
    dados(i).height = holder.ans(2,:);
    
    i = i+1;
end

%% Corte geral
for c=1:size(dados,2)
    
    %% Plot do resultado
    temp_fig = figure();
    plot(dados(c).tempo,dados(c).height);
    
    CORTE_START = input("Corte inicial - t= ");
    CORTE_END = input("Corte final - t= ");
    CORTE_START_COMMAND = input("Comando inicio - t= ");
    CORTE_K = input("Comando k - h= ");
    
    shearch = dados(c).tempo - CORTE_START;
    [~,CORTE_START] = min(abs(shearch));
    
    shearch = dados(c).tempo - CORTE_END;
    [~,CORTE_END] = min(abs(shearch));
    
    shearch = dados(c).tempo - CORTE_START_COMMAND;
    [~,CORTE_START_COMMAND] = min(abs(shearch));
    
    
    dados(c).tempo      = transpose(dados(c).tempo(CORTE_START:CORTE_END));
    dados(c).height     = transpose(dados(c).height(CORTE_START:CORTE_END) - 0.75);
    dados(c).commands   = transpose([zeros(1,CORTE_START_COMMAND-CORTE_START),...
                           ones(1,CORTE_END-CORTE_START_COMMAND+1)*CORTE_K]);
    figure();
    plot(dados(c).tempo,dados(c).height);
    hold on;
    plot(dados(c).tempo,dados(c).commands);
    title(dados(c).nome);
    
    close(temp_fig);
end

systemIdentification;