clear all;
close all;

%NOTA: para este código funcionar:
%
%      Nome dos ficheiros: 
%
%      Devem, preferencialmente, apenas diferir no sufixo co os comandos a
%      terem *_commands.mat
%      

path(pathdef);
path = [pwd '\Dados\Theta\'];

i=1;
previous = 'lol';

%% Carregar dados
files = dir(strcat(path,'*.mat'));
for file = files'
    dados(i).nome   = file.name;
    dados(i).KH     = erase(file.name,"_commands.mat");
    dados(i).KH     = erase(dados(i).KH,".mat");
    dados(i).KH     = strrep(dados(i).KH,'_',' ');
    holder          = load(strcat(path,file.name));
    if(size(holder.ans,1) == 2)             % Caso venha do vetor de comando
        dados(i).tempo  = holder.ans(1,:);
        dados(i).theta  = holder.ans(2,:);
    else                                    % Caso venha do vetor de estados
        dados(i).theta_exp  = holder.ans(3,:);
        i = i - 1;
    end
    i = i+1;
end

%% Analise geral
for c=1:size(dados,2)
    
    %% Analise do comando
    
    % Estabelece os intervalos onde os comandos estao a funcionar
    command_v  = [0,0];
    command_n   = 0;
    for i=2:(size(dados(c).tempo,2)-1)
        command_v(2) = command_on(dados(c).theta(i-1),dados(c).theta(i),dados(c).theta(i+1));
        % Caso o comando tenha começado
        if(command_v(2) == 1 && command_v(1) == 0)
            command_n = command_n +1;
            dados(c).command(command_n).start_ID = i;
        end
        % Caso o comando tenha acabado
        if(command_v(2) == 0 && command_v(1) == 1)
            dados(c).command(command_n).end_ID = i;
        end
        command_v(1) = command_v(2);
    end
    % Estabelece os intervalos onde os resultados vao ser analisados
    slack_start = 0.05;
    slack_end   = 0.05;
    for i=1:command_n
        start_c = dados(c).command(i).start_ID;
        end_c   = dados(c).command(i).end_ID;
        dados(c).analysis(i).start_ID   = start_c + ceil((end_c - start_c)*slack_start);
        dados(c).analysis(i).end_ID     = end_c   - ceil((end_c - start_c)*slack_start);
    end
    
    % Estabelece as amplitudes e frequencias dos comandos
    for i=1:command_n
        maximum_m = 0;
        maximum_n = 0;
        zeros_n   = 0;
        
        zeros_t     = zeros(1,100);
        maximum_ID  = zeros(1,100);
        % Apanha os zeros e maximos/minimos
        for j=(dados(c).analysis(i).start_ID+1):(dados(c).analysis(i).end_ID-1)
            x1 = [dados(c).tempo(j-1),dados(c).theta(j-1)];
            x2 = [dados(c).tempo(j),dados(c).theta(j)];
            x3 = [dados(c).tempo(j+1),dados(c).theta(j+1)];
            
            t = zerocross_detection(x1,x2);
            if(t~=0)
                zeros_n         = zeros_n +1;
                zeros_t(zeros_n)  = t;
            end
            if(local_maximum(abs(x1),abs(x2),abs(x3))~=0)
                maximum_n               = maximum_n +1;
                maximum_ID(maximum_n)   = j;
            end
        end
        % Calcula a frequencia
        T = 2*(zeros_t(zeros_n) - zeros_t(1))/(zeros_n-1);
        frequency = 1/T;
        % Calcula a amplitude
        amplitude = mean(abs(dados(c).theta(maximum_ID(1:maximum_n))));
        
        dados(c).command(i).zeros_t     = zeros_t(1:zeros_n);
        dados(c).command(i).maximum_ID  = maximum_ID(1:maximum_n);
        dados(c).command(i).frequency   = frequency*2*pi;
        dados(c).command(i).amplitude   = amplitude;
    end
    
    %% Analise experimental
    
    % Estabelece as amplitudes e frequencias das experiencias
    for i=1:command_n
        maximum_m       = 0;
        maximum_n       = 0;
        zeros_n         = [0,0];
        
        zeros_t         = zeros(1,100);
        maximum_ID      = zeros(1,100);
        maximum_m_ID    = zeros(1,100);
        % Apanha os zeros e maximos/minimos
        for j=(dados(c).analysis(i).start_ID+1):(dados(c).analysis(i).end_ID-1)
            x1 = [dados(c).tempo(j-1)   ,dados(c).theta_exp(j-1)];
            x2 = [dados(c).tempo(j)     ,dados(c).theta_exp(j)];
            x3 = [dados(c).tempo(j+1)   ,dados(c).theta_exp(j+1)];
            
            t = zerocross_detection(x1,x2);
            if(t~=0 && zeros_n(2)<=100)
                zeros_n(2)          = zeros_n(2) +1;
                zeros_t(zeros_n(2)) = t;
            end
            % extremos locais (todos)
            if(zeros_n(1)~=0)
                if(local_maximum(abs(x1),abs(x2),abs(x3))~=0 && maximum_m<=100)
                    maximum_m               = maximum_m +1;
                    maximum_m_ID(maximum_m) = j;
                end
                % extremo local maximo
                if(zeros_n(1)~=zeros_n(2) && maximum_n<=100 && maximum_m > 0)
                    maximum_n               = maximum_n +1;
                    position = abs_maximum(abs(dados(c).theta_exp(maximum_m_ID(1:maximum_m))));
                    maximum_ID(maximum_n)   = maximum_m_ID(position);

                    maximum_m       = 0;
                end
            end
            zeros_n(1) = zeros_n(2);
        end
        % Calcula a frequencia
        T = 2*(zeros_t(zeros_n) - zeros_t(1))/(zeros_n-1);
        frequency = 1/T;
        % Calcula a amplitude e ganho
        amplitude = mean(abs(dados(c).theta_exp(maximum_ID(1:maximum_n))));
        gain      = 20*log(amplitude/dados(c).command(i).amplitude)/log(10);
        
        dados(c).experiment(i).zeros_t      = zeros_t(1:zeros_n(1));
        dados(c).experiment(i).maximum_ID   = maximum_ID(1:maximum_n);
        dados(c).experiment(i).frequency    = frequency*2*pi;
        dados(c).experiment(i).amplitude    = amplitude;
        dados(c).experiment(i).gain         = gain;
    end
    
end

%% Estimacao da funcao de transferencia

% Resultados experimentais
plus = 0;
for c=1:size(dados,2)
    for i=1:size(dados(c).command,2)
        Transfer_function(plus + i,1) = dados(c).command(i).frequency;
        Transfer_function(plus + i,2)  = dados(c).experiment(i).gain;
        Transfer_function(plus + i,3)  = dados(c).command(i).amplitude;
        Transfer_function(plus + i,4)  = dados(c).experiment(i).amplitude;
        
    end
    plus = plus + i;
end
[~,I] = sort(Transfer_function);
Transfer_function = Transfer_function(I(:,1),[1,2,3,4]);
s = tf('s');

%% Plot dos resultados
for coisa = dados
    
    plot_command    = 1;
    plot_experiment = 1;
    figure();
    hold on;
    if(plot_command == 1)
        plot(coisa.tempo,coisa.theta);
    end
    hold on;
    if(plot_experiment == 1)
        plot(coisa.tempo,coisa.theta_exp);
    end
    % Plot dos inicios e fins de analise
    for i=1:size(coisa.command,2)
        plot([coisa.tempo(coisa.analysis(i).start_ID),coisa.tempo(coisa.analysis(i).start_ID)],...
             0.3*[-1,1],...
             '--','Color','black');
        plot([coisa.tempo(coisa.analysis(i).end_ID),coisa.tempo(coisa.analysis(i).end_ID)],...
             0.3*[-1,1],...
             '--','Color','black');
    end
    % Plot pontos comando
    if(0 == 1)
        % Plot dos zeros de comando
        for i=1:size(coisa.command,2)
            for j=1:size(coisa.command(i).zeros_t,2)
                plot(coisa.command(i).zeros_t(j),...
                     0,...
                     'ro');
            end
        end
        % Plot dos maximos de comando
        for i=1:size(coisa.command,2)
            for j=1:size(coisa.command(i).maximum_ID,2)
                plot(coisa.tempo(coisa.command(i).maximum_ID(j)),...
                     coisa.theta(coisa.command(i).maximum_ID(j)),...
                     'go');
            end
        end
    end
    
    % Plot pontos experiencia
    if(plot_experiment == 1)
        % Plot dos zeros de experiencia
        for i=1:size(coisa.experiment,2)
            for j=1:size(coisa.experiment(i).zeros_t,2)
                plot(coisa.experiment(i).zeros_t(j),...
                     0,...
                     'rx');
            end
        end
        % Plot dos maximos de experiencia
        for i=1:size(coisa.experiment,2)
            for j=1:size(coisa.experiment(i).maximum_ID,2)
                plot(coisa.tempo(coisa.experiment(i).maximum_ID(j)),...
                     coisa.theta_exp(coisa.experiment(i).maximum_ID(j)),...
                     'gx');
            end
        end
    end

    title(coisa.KH);
    legend('comando', 'experimental', 'start', 'end', 'Location', 'southeast');
    
end

%% Bode plot - Segunda ordem
fprintf('Resultados para o fitting para uma FT de <strong>2a ordem sem zero</strong>\n\n');

beta_0a = [5;10;0.3];
% Estimacao com os dados
[beta_a, ~, ~, ~, error, ~] = nlinfit(Transfer_function(:,1),...
                                    10.^(Transfer_function(:,2)/20),...
                                    @tf_pitch_a,...
                                    beta_0a...
                                    );
beta_a(2) = abs(beta_a(2));
fprintf(' K_p = %f \n omega_n = %f rad/s \n csi = %f\n', beta_a(1), beta_a(2), beta_a(3));
fprintf(' MSE - %.3f 10^-2 \n\n', error*100);

tf_a = beta_a(1)/(s^2 + 2*beta_a(2)*beta_a(3)*s + beta_a(2)^2);

opts = bodeoptions;
opts.XLim = 10.^[-1,2];

figure();
hold on;
semilogx(Transfer_function(:,1),Transfer_function(:,2),'*','MarkerEdgeColor',"#FF6100");
bodeplot(tf_a,opts);

%% Bode plot - Segunda ordem + Zero

fprintf('Resultados para o fitting para uma FT de <strong>2a ordem com zero</strong>\n\n');

beta_0b = [5;0.5;10;0.3];
% Estimacao com os dados
[beta_b, ~, ~, ~, error, ~] = nlinfit(Transfer_function(:,1),...
                                    10.^(Transfer_function(:,2)/20),...
                                    @tf_pitch_b,...
                                    beta_0b...
                                    );
beta_b(3) = abs(beta_b(3));
fprintf(' K_p = %f \n a = %f \n omega_n = %f rad/s \n csi = %f \n', ...
        beta_b(1), beta_b(2), beta_b(3), beta_b(4));
fprintf(' MSE - %.3f 10^-2 \n\n', error*100);

tf_b = beta_b(1)*(s + beta_b(2))/(s^2 + 2*beta_b(3)*beta_b(4)*s + beta_b(3)^2);

figure();
hold on;
semilogx(Transfer_function(:,1),Transfer_function(:,2),'*','MarkerEdgeColor',"#FF6100");
bodeplot(tf_b,opts);

%% Bode plot - Terceira ordem + Zero

fprintf('Resultados para o fitting para uma FT de <strong>3a ordem com zero</strong>\n\n');
beta_0c = [5;0.5;3;0.3;3];
% Estimacao com os dados
[beta_c, ~, ~, ~, error, ~] = nlinfit(Transfer_function(:,1),...
                                    10.^(Transfer_function(:,2)/20),...
                                    @tf_pitch_c,...
                                    beta_0c...
                                    );
                                
beta_c(3) = abs(beta_c(3));
fprintf(' K_p = %f \n a = %f \n omega_n = %f rad/s \n csi = %f \n b = %f \n', ...
        beta_c(1), beta_c(2), beta_c(3), beta_c(4), beta_c(5));
fprintf(' MSE - %.3f 10^-2 \n\n', error*100);

tf_c = beta_c(1)*(s + beta_c(2))/(s^2 + 2*beta_c(3)*beta_c(4)*s + beta_c(3)^2)/(s + beta_c(5));

figure();
hold on;
semilogx(Transfer_function(:,1),Transfer_function(:,2),'*','MarkerEdgeColor',"#267F00");
bodeplot(tf_c,opts);

figure();
hold on;
semilogx(Transfer_function(:,1),Transfer_function(:,2),'*','MarkerEdgeColor',"#267F00");
bodeplot(tf_a,opts);
bodeplot(tf_b,opts);


%% Funcoes
function t = zerocross_detection(x1, x2)

    if(x1(2)*x2(2) > 0)
        t = 0;
    else
        if(abs(x1(2) - x2(2)) < 10^(-6))
            t = 0;
        else
            t = x1(1) - x1(2)*((x2(1) - x1(1))/(x2(2) - x1(2)));
        end
    end
end

function i = command_on(x1, x2, x3)

    if(abs(x1 + x2 + x3) > 10^(-6))
        i = 1;
    else
        i = 0;
    end

end

function i = local_maximum(x1, x2, x3)

    if( (x1(2) <= x2(2)) && (x3(2) < x2(2)) )
        i = 1;
    else
        i = 0;
    end

end

function x = abs_maximum(list_x)

    [~, I] = max(list_x,[],2);
    
    x = I;

end

function y = tf_pitch_a(beta, x)

    % beta(1) = kp
    % beta(2) = omega_n
    % beta(3) = csi
    
    y = beta(1)./( sqrt( (beta(2)^2-x.^2).^2 + (2*beta(2)*beta(3).*x).^2 ) );
    
end
function y = tf_pitch_b(beta, x)

    % beta(1) = kp
    % beta(2) = a
    % beta(3) = omega_n
    % beta(4) = csi
    
    y = ( beta(1).*sqrt(x.^2 + beta(2)^2) )./( sqrt( (beta(3)^2-x.^2).^2 + (2*beta(3)*beta(4).*x).^2 ) );
    
end
function y = tf_pitch_c(beta, x)

    % beta(1) = kp
    % beta(2) = a
    % beta(3) = omega_n
    % beta(4) = csi
    % beta(5) = b
    
    %y = abs(beta(2)*(1i*x + beta(1)) ./ ( (1i.*x + beta(3)).*(1i.*x + beta(4)) ));
    y = ( beta(1).*sqrt(x.^2 + beta(2)^2) )./( sqrt( (beta(3)^2-x.^2).^2 + (2*beta(3)*beta(4).*x).^2 ) )./sqrt(x.^2 + beta(5)^2);
    
end