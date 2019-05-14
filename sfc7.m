function sfc7(t0,tf,dn,d0,armo,a,b,hFig)
% t0 el valor inicial para calcular la serie
% tf el valor final donde calcular la serie
% dn funci�n de la f�rmula de los dn
% f funci�n original
% armo n�mero de armonicos a utilizar en la gr�fica
% a, b intevalo para realizar la grafica de la serie

w0=2*pi/(tf-t0);

sf=d0;
t=a:0.0001:b;

for n=1:armo
    sf=sf+dn*exp(w0*-n*t*j)+dn*exp(w0*n*t*j);
end

set(hFig, 'Position', [0 0 900 900])
subplot(3,2,1)
plot(t,sf,'LineWidth',2)
grid on
legend('Serie de Fourier','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)
nn=-armo:armo;
sf=d0;
t1=t0:0.0001:tf;

for n=1:armo
    sf=sf+dn*exp(w0*-n*t1*j)+dn*exp(w0*n*t1*j);
end

absdn=zeros(1,length(nn));
cont=1;
for i =-armo:armo
    if i==0
        absdn(cont)=d0;
    end
    
    absdn(cont)=dn;
    cont=cont+1;
end

subplot(3,2,2)
stem(w0*nn,abs(absdn),'LineWidth',2)
title('Espectro de magnitud D_n ','FontWeight','bold','FontSize',16)
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on

subplot(3,2,3) % % 
stem(w0*nn,angle(absdn),'LineWidth',2) % % 
title('Espectro de fase, \angle de D_n ','FontWeight','bold','FontSize',16) % % 
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on

end