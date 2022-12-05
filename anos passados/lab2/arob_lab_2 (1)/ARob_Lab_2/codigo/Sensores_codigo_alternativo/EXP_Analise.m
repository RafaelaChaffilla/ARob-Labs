clc
clear all
close all

g = 9.8;
k_a = 0.009622374747045;
k_gir = 0.015688537908311;
k_h_raw = 1.145/(1.214e-42);
EXP = 'D'; 
switch EXP
    case 'A' %EXP A
        load('sensor_A.mat');
        load('states_A.mat');
 
    case 'B' %EXP B
        load('sensor_B.mat');
        load('states_B.mat');
    case 'C' %EXP C
        load('sensor_C.mat');
        load('states_C.mat');
    case 'D' %EXP D
        load('sensor_D.mat');
        load('states_D.mat');
        
end


t = sensor.time;
a_x=k_a*sensor.signals.values(:,1);
a_y=k_a*sensor.signals.values(:,2);
a_z=k_a*sensor.signals.values(:,3);


ver = 1;
if ver == 1
    figure (1)
    plot(t,a_x,'r');
    hold on
    plot(t,a_y,'g');
    hold on
    plot(t,a_z,'b');
    xlabel('t[s]')
    ylabel('a[N]')
    legend('a_x','a_y','a_z')
    title('Aceleração')
end

g_data_x=180/pi*k_gir*sensor.signals.values(:,4);
g_data_y=180/pi*k_gir*sensor.signals.values(:,5);
g_data_z=180/pi*k_gir*sensor.signals.values(:,6);



g_int_x=zeros(1,length(g_data_x));
g_int_y=zeros(1,length(g_data_x));
g_int_z=zeros(1,length(g_data_x));
for i = 2:length(g_data_x)
    if g_data_z(i) > 20
        g_data_z(i) = 0;
    end
    g_int_x(i)=(g_data_x(i)+g_data_x(i-1))/2*(t(i)-t(i-1))+ g_int_x(i-1);
    g_int_y(i)=(g_data_y(i)+g_data_y(i-1))/2*(t(i)-t(i-1))+ g_int_y(i-1);
    g_int_z(i)=(g_data_z(i)+g_data_z(i-1))/2*(t(i)-t(i-1))+ g_int_z(i-1);
end


rol=rad2deg(states.signals.values(:,1));
pi=rad2deg(states.signals.values(:,2));
yaw=rad2deg(states.signals.values(:,3));
t_states = states.time;

ver = 1;
if ver == 1
    figure (2)
    plot(t,g_data_x,'r')
    hold on
    plot(t,g_data_y,'g')
    hold on
    plot(t,g_data_z,'b')    
    xlabel('t[s]')
    ylabel('velocidade angular[º/s]')
    legend('p','q','r')
    title('Velocidade angular medida nos Giroscópios')
    %ylim([-1.5 1.5])
    
    
    
    figure (3)
    plot(t,g_int_x,'r')
    hold on
    plot(t_states,rol,'k')
    xlabel('t[s]')
    ylabel('Roll[º]')
    legend('Roll Sensor','Roll Real')
    title('Roll') 
    
    figure (4)
    plot(t,g_int_y,'r')
    hold on
    plot(t_states,pi,'k')
    xlabel('t[s]')
    ylabel('Pitch[º]')
    legend('Pitch Sensor','Pitch Real')
    title('Pitch') 

    %ylim([-1.5 1.5])
    

    figure (5)
    plot(t,g_int_z,'r')
    hold on
    plot(t_states,yaw,'k')
    xlabel('t[s]')
    ylabel('Yaw[º]')
    legend('Yaw Sensor','Yaw Real')
    title('Yaw') 
    

end

h=1/0.75*states.signals.values(:,9);

ver = 1;


ver = 1;
if ver == 1
    figure (7)
    plot(t,g_int_x,'r');
    hold on
    plot(t,g_int_y,'g');
    hold on
    plot(t,g_int_z,'b');
    hold on
    xlabel('t[s]')
    ylabel('Posição Angular[º]')
    legend('rho','pitch','yaw')
    title('Angulos de Euler')
end




meanvec=[mean(a_x(2000:5000)),mean(a_y(2000:5000)),mean(a_z(2000:5000)),mean(g_data_x(2000:5000)),mean(g_data_y(2000:5000)),mean(g_data_z(2000:5000))];

covvec=[cov(a_x(2000:5000)),cov(a_y(2000:5000)),cov(a_z(2000:5000)),cov(g_data_x(2000:5000)),cov(g_data_y(2000:5000)),cov(g_data_z(2000:5000))];
        

plot(t,sensor.signals.values(:,1),'r')
hold on
plot(t,sensor.signals.values(:,2),'g')
hold on
plot(t,sensor.signals.values(:,3),'b')
    xlabel('t[s]')
    ylabel('a')
    legend('asen_x','asen_y','asen_z')
    title('Dados dos Acelerômetros')
    
    
    
    
    subplot(1,2,1)
    plot(t,1/k_gir*g_int_x,'r')
    xlabel('t[s]')
    ylabel('rho')
    title('Giroscópio x_B int')
    subplot(1,2,2)
    plot(t_states,rol,'b')
    xlabel('t[s]')
    ylabel('rho[º]')
    title('Rolamento')
    

    
    
    
    plot(t,g_data_x,'r')
    hold on
    plot(t,g_data_y,'g')
    hold on
    plot(t,g_data_z,'b')
    hold on
    %ylim([-50 50])
    xlabel('t[s]')
    ylabel('Velocidade Angular[º/s]')
    legend('Giroscópio X','Giroscópio Y','Giroscópio Z')
    title('Medições dos giroscópios - Exp B')
    
    
    plot(t,a_x,'r')
hold on
plot(t,a_y,'g')
hold on
plot(t,a_z,'b')
    xlabel('t[s]')
    ylabel('a[m/s^2]')
    legend('a_x','a_y','a_z')
    title('Acelerômetros - EXP C')
    
    


    
        plot(t,g_data_x,'r')
    hold on
    plot(t,g_data_y,'g')
    hold on
    plot(t,g_data_z,'b')
    hold on
    %ylim([-50 50])
    xlabel('t[s]')
    ylabel('Velocidade Angular[º/s]')
    legend('Giroscópio X','Giroscópio Y','Giroscópio Z')
    title('Medições dos giroscópios - Exp A')
    
    
    figure(10)
    subplot(1,2,1)
    plot(t,rad2deg(1/g*a_x-0.012),'r')
    hold on
    plot(t_states,pi,'k')
    %ylim([-50 50])
    xlabel('t[s]')
    ylabel('Pitch[º]')
    legend('Pitch Inclinómetro','Pitch Estado')
    title('Inclinómetro Pitch - Exp D')
   
    subplot(1,2,2)
    plot(t, rad2deg(-1/g*a_y-0.032),'b')
    hold on
    plot(t_states,rol,'k')
    %ylim([-50 50])
    xlabel('t[s]')
    ylabel('Rol[º]')
    legend('Rol Inclinómetro','Rol Estado')
    title('Inclinómetro Rol - Exp D')
    
    
    
    
close all
    
   
    
    