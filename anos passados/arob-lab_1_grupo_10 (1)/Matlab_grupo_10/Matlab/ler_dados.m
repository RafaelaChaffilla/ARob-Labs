%Lê os dados guardados nos ficheiros de texto

fileID = fopen('dados4.txt', 'r');

A = fscanf(fileID, '%f   %f', [2 Inf]);

A = A';

%Vetor C só com os dados relevantes para o systemIdentification - decidido manualmente.
%De seguida é removido o valor do primeiro step.
C = A(232:663, 2);
C = C - 1;

plot(A(:,1),A(:,2))

%Vetor de input para o System Identification
for i = 1:432
    B(i) = 0.2;
end
 
 B = B';

%plot das resposta da picada - não utilizado
%figure;
% plot(tout, ref.signals.values, tout, simout.signals.values, A(:,1), A(:,3));
% title('Resposta a uma referência')
% xlabel('Tempo (s)');
% ylabel('Ângulo de picada (rad)');
% legend('Referência', 'Simulação', 'Experimental');
