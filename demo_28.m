Z_R = 30-40i;
Z_0 = 50;
lamda = 1;
Gama_R = (Z_R-Z_0)/(Z_R+Z_0);
b = 2*abs(Gama_R)/(1-abs(Gama_R)^2)^(1/2);
b1 = b;
b2 = -b;
ds1 = lamda/(4*pi)*(angle(Gama_R)-pi/2-atan(b1/2));
n = 1;
 while ds1<0 || ds1>=lamda/2
        ds1 = lamda/(4*pi)*(angle(Gama_R)-pi/2-atan(b1/2)+2*n*pi);
        n = n+1;
 end
 ls1 = lamda/(2*pi)*atan(-1/b1)+lamda/2;
 n = 1;
 ds2 = lamda/(4*pi)*(angle(Gama_R)+pi/2-atan(b2/2));
 while ds2<0 || ds2>=lamda/2
        ds2 = lamda/(4*pi)*(angle(Gama_R)+pi/2-atan(b2/2)+2*n*pi)
        n = n+1;
 end
 ls2 = lamda/(2*pi)*atan(-1/b2);

 vp = 3*10^8;
 T = lamda/vp;
 omega = 2*pi*vp/lamda;
 beta = 2*pi/lamda;
 A = 5;
 V_0 = A;
Gama2 = 0;
Gama1 = Gama_R;
Gama3 = -1;
z2 = -1:0.01:0;
z1 = -ds1:ds1/100:0;
z3 = -ls1:ls1/100:0;
V2_plus = V_0;
V1_plus = V_0 / (exp(i*beta*ds1) .* (1+Gama1*exp(-i*2*beta*ds1)));
V3_plus = V_0 / (exp(i*beta*ls1) .* (1+Gama3*exp(-i*2*beta*ls1)));
ph_V2 = V2_plus * exp(-i*beta*z2);
ph_V1 = V1_plus * exp(-i*beta*z1) .* (1+Gama1*exp(i*2*beta*z1));
ph_V3 = V3_plus * exp(-i*beta*z3) .* (1+Gama3*exp(i*2*beta*z3));
ph_I2 = V2_plus/Z_0 * exp(-i*beta*z2);
ph_I1 = V1_plus/Z_0 * exp(-i*beta*z1) .* (1-Gama1*exp(i*2*beta*z1));
ph_I3 = V3_plus/Z_0 * exp(-i*beta*z3) .* (1-Gama3*exp(i*2*beta*z3));

z2_prime = z2 - ds1;
ph_V12 = cat(2, ph_V2, ph_V1);
ph_I12 = cat(2, ph_I2, ph_I1);
z12 = cat(2, z2_prime, z1);
V2 = real(ph_V2);
V1 = real(ph_V1);
V3 = real(ph_V3);
I2 = real(ph_I2);
I1 = real(ph_I1);
I3 = real(ph_I3);
V12 = real(ph_V12);
I12 = real(ph_I12);

figure
h = plot(z3, V3, 'EraseMode', 'xor');
xlabel('z3 (m)');
ylabel('V3 (Volt)');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'position',[-ls1 8.5])
axis([-ls1 0 -8 8]);			% 設定圖軸的範圍
grid on					% 畫出格線
tic
for j = 0:1000
    t = j/1000 * T;
    y = real(ph_V3 * exp(i*omega*t));
    set(h, 'ydata', y);		% 設定新的 y 座標
	drawnow				% 立即作圖
end
toc