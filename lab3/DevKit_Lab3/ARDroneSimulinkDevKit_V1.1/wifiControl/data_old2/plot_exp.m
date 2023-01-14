close all;
clear all;

%fixme - colocar dados certos
Alt_exp.desired = load('desired_data_3.mat');
Alt_exp.results = load('exp_data_3.mat');

% Plot altitude tracking
fig1 = figure();
hold on;
plot(Alt_exp.desired);
for c = 1:size(sim_results,2)
    plot(sim_results(c).tout, sim_results(c).h_sim.signals.values);
end
xlim([0 30]);
legend('pedida', sprintf('k = %0.1f', k_w_vector(1)), ...
    sprintf('k_w = %0.1f', k_w_vector(2)), ...
    sprintf('k_w = %0.1f', k_w_vector(3)), ...
    sprintf('k_w = %0.1f', k_w_vector(4)), ...
    sprintf('k_w = %0.1f', k_w_vector(5)), ...
    sprintf('k_w = %0.1f', k_w_vector(6)), ...
    'location', 'southeast');
xlabel('Tempo [s]');
ylabel('Altitude h [m]');
saveas(fig1, '../imgs/Altitude_Sim', 'png');