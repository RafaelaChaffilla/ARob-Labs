clear all;
close all;

%NOTA: para este código funcionar:
%
%      Nome dos ficheiros: 
%
%      1 -  devem começar com a indicação do k e do h, e deve estar
%      exatamente da mesma forma escrito para o experimental e o teorico
%      
%      2 -  depois do k e h, dizer "_teorico", "_experimental"
%

path(pathdef);
path_1 = [pwd '\Dados\H\Raw\'];
path_2 = [pwd '\Dados\H\'];
addpath([pwd '\Functions']);

i=1;
previous = 'lol';

if isfile('H_DADOS.mat')
    load('H_DADOS.mat')
else
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
end

if isfile('estimated_tf_height.mat')
    load('estimated_tf_height.mat')
else
    systemIdentification;
end

color_new =  [  0           0.4470    0.7410
                0.8500      0.3250    0.0980
                0.9290      0.6940    0.1250
                0.4940      0.1840    0.5560
                0.4660      0.6740    0.1880
                0.3010      0.7450    0.9330
                0.6350      0.0780    0.1840];

figure();
hold on;
grid('on');
xlabel('Real Axis [s^{-1}]');
ylabel('Imaginary Axis [s^{-1}]');
axis([-3 1 -3 3]);
for i=1:size(tf_master,1)
    
    z_1 =   -tf_master(i).Denominator(2)/2 + ...
            sqrt(tf_master(i).Denominator(2)^2 - 4*tf_master(i).Denominator(3))/2;
    z_2 =   -tf_master(i).Denominator(2)/2 - ...
            sqrt(tf_master(i).Denominator(2)^2 - 4*tf_master(i).Denominator(3))/2;
    
    if(i == 1)
        j = 0.5;
    else
        j = i-1;
    end
    name = ['K_p = ' num2str(j)];
    
    h(2*i-1) = plot(real(z_1),imag(z_1),'*','Color',color_new(i,:),'DisplayName',name);
    h(2*i-1) = plot(real(z_2),imag(z_2),'*','Color',color_new(i,:),'DisplayName',name);

end


legend([h(1),h(3),h(5),h(7)]);
