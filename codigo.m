A=3;
f=@(t) (2*A*t.*(t>=-1/2 & t<=1/2))+(2*A*(1-t).*(t>1/2 & t<=3/2));
t0=-1/2;%valor inicial para hacer la serie
tf=3/2;%valor final para hacer la serie
t=t0:0.0001:tf;
sumterms = zeros(16, length(t)); 
sumterms(1,:) = 1/2;
for n = 1:size(sumterms,1)-1
sumterms(n+1,:) = (2/(pi*n)*sin(pi*n/2))*cos(n*t);
end
x_N = cumsum (sumterms);
figure(14); 
clf; 
ind = 0;
for N = [0,1:2:size(sumterms, 1)-1]
ind = ind+1;
subplot (3,3,ind);
plot (t,x_N(N+1), 'k',t,f(t), 'k--') 
axis ([-2*pi 2*pi -0.2 1.2]);
xlabel ('t'); 
ylabel (['x_{',num2str(N),'} (t)']);
end