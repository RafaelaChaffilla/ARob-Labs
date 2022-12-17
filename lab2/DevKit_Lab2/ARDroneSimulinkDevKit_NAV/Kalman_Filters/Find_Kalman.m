%% Find the best Kalman Coeficients

%% Definitions
clear all;

MAX_GEN                     = 20;
MAX_GEN_WITHOUT_IMPROVEMENT = 10;
SIZE_POP                    = 20;
BEST_BOYS                   = 4;
SONS                        = (SIZE_POP-BEST_BOYS)/BEST_BOYS;
conter_MGWI = 0;

%% Loads
files = ['data_for_tunning\D.mat';...
         'data_for_tunning\E.mat'];
     
for c=1:size(files,1)
    Temp = load(files(c,:));
    Temp = Temp.ans;
    DATA(c).accs_meas = [Temp(1,:); Temp(2:4,:)];
    DATA(c).gyro_meas = [Temp(1,:); Temp(5:7,:)];
    DATA(c).angs_real = [Temp(1,:); Temp(8:9,:)];
end

%% Initial Values

population(1).Q_gyro = 1;
population(1).Q_bias = 10^(-6);
population(1).R_accs = 1;

history(1:MAX_GEN+1) = population(1);
%% Define and Evaluate Pops

% Gen 0
for c=2:SIZE_POP
    population(c) = scramble(population(1));
    fprintf("%3d Kalmans made.\n",c);
end

% Gen N
for c=0:MAX_GEN
    % Scrambles the eggs (if it isn't Gen 0)
    if(c~=0)
        for d=1:BEST_BOYS
            for e=BEST_BOYS+SONS*(d-1)+1:BEST_BOYS+SONS*d
                population(e) = scramble(population(d));
            end
        end
    end
    % Rates them
    for d=1:size(population,2)
        results(d) = evaluation(population(d), DATA);
        fprintf("%3d Kalmans rated - %f.\n",d,results(d));
    end
    % Sorts
    [~,I] = sort(results, 'ascend');
    population = population(I);    
    % Verifies Conditions
    if I(1) == 1
        conter_MGWI = conter_MGWI +1;
    else
        conter_MGWI = 0;
    end
    if conter_MGWI >= MAX_GEN_WITHOUT_IMPROVEMENT
       break; 
    end
    history(c+2) = population(1);
    fprintf("%d gen. Times without change %d.\n",c,conter_MGWI);
end
for c = 1:size(history,2)
    History.Q_gyro(c) = history(c).Q_gyro;
    History.Q_bias(c) = history(c).Q_bias;
    History.R_accs(c) = history(c).R_accs;
end
figure();
scatter(History.Q_bias,History.R_accs);
figure();
comet(History.Q_bias,History.R_accs);
%% Final Results
Kalman = population(1);
Q           = [ones(3,1)*Kalman.Q_gyro;...
               ones(3,1)*Kalman.Q_bias];
R           =  ones(3,1)*Kalman.R_accs;
Kalman_OP   = Setup_Kalman_OP(Q, R);
for c = 1:size(DATA,2)
    SIM.angs_real = DATA(c).angs_real';
    SIM.gyro_meas = DATA(c).gyro_meas';
    SIM.accs_meas = DATA(c).accs_meas';
    
    finalTime   = SIM.angs_real(end,1);
    sampleTime  = SIM.angs_real(2,1)-SIM.angs_real(1,1);
    
    options     = simset('SrcWorkspace','current');
    rrr(c)  = sim('Kalman_Tuning',[],options);

end

%% Functions

% Scrambler
function son = scramble(dad)
% Random number from 0 - 0.1
luck = rand(2,2);
son.Q_gyro = dad.Q_gyro;
son.Q_bias = dad.Q_bias*(1 + sign(luck(1,2) - 0.5)*luck(1,1)*0.1);
son.R_accs = dad.R_accs*(1 + sign(luck(2,2) - 0.5)*luck(2,1)*0.1);
end
% Evaluation

function score = evaluation(Kalman, DATA)
% Prepares the Kalman Filter
Q           = [ones(3,1)*Kalman.Q_gyro;...
               ones(3,1)*Kalman.Q_bias];
R           =  ones(3,1)*Kalman.R_accs;
Kalman_OP   = Setup_Kalman_OP(Q, R);

score = 0;
for c = 1:size(DATA,2)
    SIM.angs_real = DATA(c).angs_real';
    SIM.gyro_meas = DATA(c).gyro_meas';
    SIM.accs_meas = DATA(c).accs_meas';
    
    finalTime   = SIM.angs_real(end,1);
    sampleTime  = SIM.angs_real(2,1)-SIM.angs_real(1,1);
    
    options = simset('SrcWorkspace','current');
    results = sim('Kalman_Tuning',[],options);
    
    score = merit(results.cost(:,2));
end
end

% function of merit
function m = merit(diference_vector)

m = 0;
for c=1:size(diference_vector,1)
    m = m + cost(diference_vector(c));
end

m = m/size(diference_vector,1);

end

% Polinomial of the cost
function c = cost(diference)
c = diference^2;
end