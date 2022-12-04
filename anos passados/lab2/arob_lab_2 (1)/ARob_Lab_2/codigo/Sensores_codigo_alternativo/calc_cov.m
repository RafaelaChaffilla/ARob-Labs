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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        

%         
%         load('sensor_A.mat')
%         sensor_K= zeros(size(sensor.signals.values));
%         sensor_K(:,1)=k_a*sensor.signals.values(:,1);
%         sensor_K(:,2)=k_a*sensor.signals.values(:,2);
%         sensor_K(:,3)=k_a*sensor.signals.values(:,3);
%         sensor_K(:,4)=k_gir*sensor.signals.values(:,4);
%         sensor_K(:,5)=k_gir*sensor.signals.values(:,5);
%         sensor_K(:,6)=k_gir*sensor.signals.values(:,6);
%         sensor_K(:,7)=k_h_vision*sensor.signals.values(:,7);
%         sensor_K(:,8)=k_h_vz*sensor.signals.values(:,8);
%         sensor_K(:,9)=k_h_raw*sensor.signals.values(:,9);
%         
%         
%         cov_A = cov(sensor_K);
%         
%         load('sensor_B.mat')
%         sensor_K= zeros(size(sensor.signals.values));
%         sensor_K(:,1)=k_a*sensor.signals.values(:,1);
%         sensor_K(:,2)=k_a*sensor.signals.values(:,2);
%         sensor_K(:,3)=k_a*sensor.signals.values(:,3);
%         sensor_K(:,4)=k_gir*sensor.signals.values(:,4);
%         sensor_K(:,5)=k_gir*sensor.signals.values(:,5);
%         sensor_K(:,6)=k_gir*sensor.signals.values(:,6);
%         sensor_K(:,7)=k_h_vision*sensor.signals.values(:,7);
%         sensor_K(:,8)=k_h_vz*sensor.signals.values(:,8);
%         sensor_K(:,9)=k_h_raw*sensor.signals.values(:,9);
%         
%         
%         cov_B = cov(sensor_Kc);
%         
        load('sensor_C.mat')
        sensor_K= zeros(size(sensor.signals.values));
        sensor_K(:,1)=k_a*sensor.signals.values(:,1);
        sensor_K(:,2)=k_a*sensor.signals.values(:,2);
        sensor_K(:,3)=k_a*sensor.signals.values(:,3);
        sensor_K(:,4)=k_gir*sensor.signals.values(:,4);
        sensor_K(:,5)=k_gir*sensor.signals.values(:,5);
        sensor_K(:,6)=k_gir*sensor.signals.values(:,6);
        sensor_K(:,7)=k_h_vision*sensor.signals.values(:,7);
        sensor_K(:,8)=k_h_vz*sensor.signals.values(:,8);
        sensor_K(:,9)=k_h_raw*sensor.signals.values(:,9);
        
        
        cov_C_acc = cov(sensor_K());
%         
%         
%         load('sensor_D.mat')
%         sensor_K= zeros(size(sensor.signals.values));
%         sensor_K(:,1)=k_a*sensor.signals.values(:,1);
%         sensor_K(:,2)=k_a*sensor.signals.values(:,2);
%         sensor_K(:,3)=k_a*sensor.signals.values(:,3);
%         sensor_K(:,4)=k_gir*sensor.signals.values(:,4);
%         sensor_K(:,5)=k_gir*sensor.signals.values(:,5);
%         sensor_K(:,6)=k_gir*sensor.signals.values(:,6);
%         sensor_K(:,7)=k_h_vision*sensor.signals.values(:,7);
%         sensor_K(:,8)=k_h_vz*sensor.signals.values(:,8);
%         sensor_K(:,9)=k_h_raw*sensor.signals.values(:,9);
%         cov_D = cov(sensor_K);
        
        
        load('sensor_A.mat')
        cov_A=cov(k_a*sensor.signals.values(2000:5000,:));
        load('sensor_B.mat')
        cov_B=cov(k_a*sensor.signals.values(2000:5000,:));
%         load('sensor_C.mat')
%         cov_C=cov(k_a*sensor.signals.values(2000:5000,:));
        load('sensor_D.mat')
        cov_D=cov(k_a*sensor.signals.values(2000:5000,:));