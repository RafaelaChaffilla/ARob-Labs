%Leitura dos dados dos ficheiros para a construção das imagens da resposta 
%aos diferentes steps e ganhos

[t1,h1] = readvars('dados1.txt');
[t2,h2] = readvars('dados2.txt');
[t3,h3] = readvars('dados3.txt');
[t4,h4] = readvars('dados4.txt');

% figure(1)
% plot(t1,h1,'g',simout.time,simout.signals.values,'r')
% hold on
% plot([0 5],[0 0],'b',[5 5],[0 1],'b',[5 15],[1 1],'b',[15 15],[1 2],'b',...
%     [15 40],[2 2],'b')
% title('Kp = 0,5')
% xlabel('Time [s]')
% ylabel('Height [m]')
% legend('Experimental','Simulação','Referência','Location','southeast')
% xlim([0 30])
% ylim([0 2.25])
% saveas(gcf,'ganho05.png')

% figure(2)
% plot(t2,h2,'g',simout.time,simout.signals.values,'r')
% hold on
% plot([0 5],[0 0],'b',[5 5],[0 1],'b',[5 15],[1 1],'b',[15 15],[1 1.5],'b',...
%     [15 40],[1.5 1.5],'b')
% title('Kp = 1')
% xlabel('Time [s]')
% ylabel('Height [m]')
% legend('Experimental','Simulação','Referência','Location','southeast')
% xlim([0 24])
% ylim([0 1.75])
% saveas(gcf,'ganho1.png')

% figure(3)
% plot(t3,h3,'g',simout.time,simout.signals.values,'r')
% hold on
% plot([0 5],[0 0],'b',[5 5],[0 1],'b',[5 15],[1 1],'b',[15 15],[1 1.25],'b',...
%     [15 40],[1.25 1.25],'b')
% title('Kp = 2,4')
% xlabel('Time [s]')
% ylabel('Height [m]')
% legend('Experimental','Simulação','Referência','Location','southeast')
% xlim([0 30])
% ylim([0 1.5])
% saveas(gcf,'ganho24.png')

figure(4)
plot(t4,h4,'g',simout.time,simout.signals.values,'r')
hold on
plot([0 5],[0 0],'b',[5 5],[0 1],'b',[5 15],[1 1],'b',[15 15],[1 1.2],'b',...
    [15 40],[1.2 1.2],'b')
title('Kp = 3,2')
xlabel('Time [s]')
ylabel('Height [m]')
legend('Experimental','Simulação','Referência','Location','southeast')
xlim([0 40])
ylim([0 1.5])
saveas(gcf,'ganho32.png')