%Script utilizado para calcular a função de transferência da picada
%usando a função nlinfit

w = [1 3 5 7 9 11 13 15 17 20];
a = [0.081 0.119 0.102 0.081 0.06 0.041 0.031 0.022 0.018 0.012];

a = a / 0.2;

my_fun = @(beta,w)(beta(1).*sqrt(w.^2+beta(2).^2)./sqrt((beta(3).^2-w.^2).^2+(2.*beta(3).*beta(4).*w).^2));

initials = [2 0.5 3 0.4];

new_coeffs = nlinfit(w, a, my_fun, initials);

disp(new_coeffs);

new_a = my_fun(new_coeffs, w);

plot(w, a, w, new_a);

disp(rms(new_a - a));