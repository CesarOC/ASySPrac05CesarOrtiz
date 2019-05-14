function sfc1(t0,tf,an,bn,a0,f,armo,a,b,hFig)
% t0 el valor inicial para calcular la serie
% tf el valor final donde calcular la serie
% dn función de la fórmula de los dn
% f función original
% armo número de armonicos a utilizar en la gráfica
% a, b intevalo para realizar la grafica de la serie

w0=2*pi/(tf-t0);

sf=a0;
t=a:0.0001:b;

for n=1:armo
     sf=sf+an(n)*cos(w0*n*t)+bn(n)*sin(w0*n*t);
end
set(hFig, 'Position', [0 0 900 900])
subplot(3,2,1)
plot(t,sf,'LineWidth',2)
grid on
legend('Serie de Fourier','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)

sf=a0;
t1=t0:0.0001:tf;

for n=1:armo
     sf=sf+an(n)*cos(w0*n*t1)+bn(n)*sin(w0*n*t1);
end

subplot(3,2,2)
plot(t1,f(t1),'r','LineWidth',2)
grid on
hold on
plot(t1,sf,'LineWidth',2)
legend('Función original','Serie de Fourier ','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)
nn=0:armo+1;
axis auto

subplot(3,2,4)
e=f(t1)-sf;
plot(t1,e,'LineWidth',2)
title('Error','FontWeight','bold','FontSize',16)
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on

subplot(3,2,6)
e=f(t1)-sf;
area(t1,e.^2)
legend('Energia del error','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on

absan=zeros(1,length(nn));
cont=1;
for i =0:armo
    if i==0
        absan(1)=0.504;
    else if i>0
    absan(cont)=an(i)*10;
    cont=cont+1;
        else
        end   
end

subplot(3,2,3)
stem(nn,absan,'LineWidth',2)
title('Valores a_n ','FontWeight','bold','FontSize',16)
xlabel('n','FontWeight','bold','FontSize',16)
grid on

absbn=zeros(1,length(nn));
conb=1;
for j =0:armo
    if j==0
        absbn(conb)=0;
    end
    absbn(conb)=bn(j);
    conb=conb+1;
end
subplot(3,2,5) % % 
stem(nn,absbn,'LineWidth',2) % % 
title('Valores b_n, b_n ','FontWeight','bold','FontSize',16) % % 
xlabel('n','FontWeight','bold','FontSize',16)
grid on

end