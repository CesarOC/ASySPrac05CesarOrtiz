syms t
n=(0:4);
f=exp(-t/2);
ti=0; tf=pi; t0=tf-ti; w0=2*pi/t0;
f2(t)=f*exp(-n*w0*1i*t);

Dn=(1/t0)*int(f2,t,ti,tf);
exacto=double (Dn);
%trapecio
m=15000;
n=(0:4);
b=pi; a=0; t=-100:.1:100;
w0=2*pi/(b-a);
f=@(t) exp(-t/2)*exp(-n*w0*1i*t);

h=(b-a)/m;

aprox=f(a)+f(b);
for i=1:m-1
x=a+i*h;
aprox=aprox+2*f(x);
end
aprox=(h/2)*aprox;
trapecio=aprox;
%metodo lathi 
T_0 = pi;
N_0 = 256; 
T = T_0/N_0; 
x = (0:T:T*(N_0-1))';

h = exp(-x/2);
h(1) = (exp(-pi/2) + 1)/2;
D_n = fft (h)/N_0;

Lathi=[D_n(1) D_n(2) D_n(3) D_n(4) D_n(5)];
%errores
error_trapecio=abs(exacto-trapecio);
error_LATHI=abs(exacto-Lathi);
%columna
coeficiente= {'d0';'d1';'d2';'d3';'d4'};
exacto=vec2mat(exacto,1);
trapecio=vec2mat(trapecio,1);
error_trapecio=vec2mat(error_trapecio,1);
Lathi=vec2mat(Lathi,1);
error_LATHI=vec2mat(error_LATHI,1);
T=table(exacto,trapecio,Lathi,error_trapecio,error_LATHI)
