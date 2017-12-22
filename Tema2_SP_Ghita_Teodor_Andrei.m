%Numar de ordine = 9
%Semnal Dreptunghiular
%Perioada P = 40 s
%Durata semnalelor D [sec] = 9
%Numar de coeficienti N = 50

P = 40;
D = 9;
N = 50;
W0 = 2*pi / P;   %Pulsatie unghiulara a semnalului dreptunghiular

pas_esant = P / 100;  
t = (-2 * P) : pas_esant : (2 * P);     %Pasul de esantionare

coef_serieF = zeros(1,N);   %Initializarea vectorului de coeficienti SFC
A = zeros(1,N);             %Initializarea vectorului de coeficienti SFA

x_t = square( W0 * t, dutyCicle); %Semnalul initial de tip dreptunghiular
x_t_init = @(t,k) square( W0 * t, dutyCicle).*exp( -1j * k * W0 * t); %Semnal descompus in SFC
x_t_reconstr = 0;    %Initializarea semnalului reconstruit
comp_cont = (1 / P) * integral(@(t) x_t_init(t,0),0,P);   %Coeficientii Xk din dezvoltarea in serie Fourier a semnalului

%Se calculeaza coeficientii SFC si se reconstruieste semnalul prin insumarea termenilor
for k = 1:50   %k reprezinta indicele termenilor din suma
    coef_serieF(k) = (1 / P) * integral(@(t) x_t_init(t,k - 25 ),0,P);   %Coeficientii Xk din dezvoltarea in serie Fourier a semnalului
    x_t_reconstr = x_t_reconstr + coef_serieF(k) * exp( 1j * (k - 25) * W0 * t);  %Recostructia semnalului cu ajutorul coeficientilor 
end                                                                               

figure(1);      %Reprezentarea semnalului dreptunghiular initial si reconstruit
subplot(2,1,1);
plot(t, x_t_reconstr, t, x_t,'r');
title('Semnalul dreptunghiular inital si reconstruit');
hold on;

A(1) = abs(comp_cont);   %Calculul coeficientiilor SFA care reprezinta amplitudiniile din spectru
for k = 1:N
      A(k+1) = 2 * abs(coef_serieF(k)); 
end

subplot(2,1,2);     %Reprezentarea spectrului de amplitudini al semnalului dreptunghiular
stem([0:N] * W0, A,'b'); 
title('Spectrul semnalului dreptunghiular');

%Conform teoriei seriilor Fourier, orice semnal periodic poate fi aproximat printr-o suma
%de sinusuri si cosinusuri de frecvente diferite si coeficienti diferiti, acesti coeficienti
%reprezentand spectrul de amplitudini al semnalului periodic initial. Prin
%folosirea unui numar finit de termeni(N=50) semnalul reconstruit pe baza
%spectrului se apropie de cel initial avand insa o marja de eroare

