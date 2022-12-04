%Covariancias

clc 
close all
clear all


Caso = 'C';
switch Caso
    case 'A' %Acelerometros
        
        load('sensor_A.mat')

        for i = 1 : size(sensor.signals.values,2)
            Average(i) = sum(sensor.signals.values(:,i))/length(sensor.signals.values(:,i));
        end

        t = sensor.time;

        a_mod_sen = sqrt(Average(1)^2+Average(2)^2+Average(3)^2);

        g=9.8;
        k_a = g/a_mod_sen;

        a_x=k_a*sensor.signals.values(:,1);
        a_y=k_a*sensor.signals.values(:,2);
        a_z=k_a*sensor.signals.values(:,3);

        ver = 1;
        if ver == 1
            Average
            figure (1)
            subplot(1,3,1)
            plot(t,a_x);
            subplot(1,3,2)
            plot(t,a_y);
            subplot(1,3,3)
            plot(t,a_z);
        end
    case 'B' %Gyroscopios
        load('sensor_D.mat')
        load('states_D.mat')
        
        t = sensor.time;
        g_data_x=sensor.signals.values(:,4);
        g_data_y=sensor.signals.values(:,5);
        g_data_z=sensor.signals.values(:,6);
        
%         a=900;
%         b=3000;
        g_int_x=zeros(1,length(g_data_x));
        g_int_y=zeros(1,length(g_data_x));
        g_int_z=zeros(1,length(g_data_x));
        for i = 2:length(g_data_x)
            g_int_x(i)=(g_data_x(i)+g_data_x(i-1))/2*(t(i)-t(i-1))+ g_int_x(i-1);
            g_int_y(i)=(g_data_y(i)+g_data_y(i-1))/2*(t(i)-t(i-1))+ g_int_y(i-1);
            g_int_z(i)=(g_data_z(i)+g_data_z(i-1))/2*(t(i)-t(i-1))+ g_int_z(i-1);
        end
        
           
        
        
        rol=states.signals.values(:,1);
        pi=states.signals.values(:,2);
        gui=states.signals.values(:,3);
        t_states = states.time;
        
        g_int_x_avg=sum(abs(g_int_x(1000:4000)))/length(g_int_x(1000:4000)); 
        rol_avg=sum(abs(rol(1000:4000)))/length(rol(1000:4000));
        
        k_gir=rol_avg/g_int_x_avg;
        ver = 0;
        if ver == 1
            figure (1)
            subplot(3,3,1)
            plot(t,g_int_x);
            title('g_x_int')
            subplot(3,3,2)
            plot(t,g_int_y);
            title('g_y_int')
            subplot(3,3,3)
            plot(t,g_int_z);
            title('g_y_int')
            subplot(3,3,4)
            plot(t,g_data_x);
            title('g_x')
            subplot(3,3,5)
            plot(t,g_data_y);
            title('g_y')
            subplot(3,3,6)
            plot(t,g_data_z);
            title('g_z')
            subplot(3,3,7)
            plot(t_states,rol);
            title('rol')
            subplot(3,3,8)
            plot(t_states,pi);
            title('pi')
            subplot(3,3,9)
            plot(t_states,gui);
            title('gui')
        end
        
        ver = 1;
        if ver == 1
            figure (1)
            plot(t_states,rol);
            hold on
            plot(t,k_gir*g_int_x);
            hold on
            title('rol')
            
            figure(2)
            plot(t_states,pi);
            hold on
            plot(t,k_gir*g_int_y);
            hold on
            title('pi')
            
            figure (3)
            plot(t_states,gui);
            hold on
            plot(t,k_gir*g_int_z);
            hold on
            title('gui')
        end
        
        
    case 'C' %Altura
        
        load('sensor_C.mat')
        load('states_C.mat')
        
          t = sensor.time;
        h_vision=sensor.signals.values(:,7);
        h_vz=sensor.signals.values(:,8);
        h_raw=sensor.signals.values(:,9);
        
        
        h=states.signals.values(:,9);
        
        h_vision_avg=sum(abs(h_vision(1000:4000)))/length(h_vision(1000:4000));
        h_vz_avg=sum(abs(h_vz(1000:4000)))/length(h_vz(1000:4000));
        h_raw_avg=sum(abs(h_raw(1000:4000)))/length(h_raw(1000:4000));
        h_avg=sum(abs(h(1000:4000)))/length(h(1000:4000));
        k_h_vision = 0.86/(1.2*10^(-42));
        k_h_vz = -0.86/533;
        k_h_raw = 0.86/(1.2*10^(-42));
        
        
        h_vz_int=zeros(1,length(h_vz));
        for i = 2:length(h_vz)
            h_vz_int(i)= (h_vz(i)+h_vz(i-1))/2*(t(i)-t(i-1))+ h_vz_int(i-1);
        end
        
        
        ver = 0;
        if ver == 1

            plot(t,h_vision);
            hond on
            plot(t,h_vz_int);
            hond on
            plot(t,h_raw);
            hond on
            plot(t,h);
            hond on
        end

        ver = 1;
        if ver == 1
           figure (2)
            plot(t,h,'b');
            title('h')
            hold on
            plot(t,k_h_vision*h_vision,'r');
            hold on
            plot(t,k_h_vz*h_vz_int,'k');
            hold on
            plot(t,k_h_raw*h_raw,'g');

            
            legend('h','h vision','h vz int','h raw')
        end
        ver = 1;
        if ver == 1
            figure (3)
            plot(t,h,'b');
            hold on
            plot(t,k_h_vision*h_vision,'r');
            legend('h','h vision')
            
            figure (4)
            plot(t,h,'b');
            hold on
            plot(t,k_h_vz*h_vz_int,'k');
            legend('h','h vz int')
            
            figure (5)
            plot(t,h,'b');
            hold on
            plot(t,k_h_raw*h_raw,'g');
            title('h raw')
            legend('h','h raw')
           
        end
end



