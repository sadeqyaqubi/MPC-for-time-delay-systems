function y = syseq (x1_k, x2_k, x3_k, x4_k,x1_k1, x2_k1, x3_k1, x4_k1, x2_kd, x2_kd1, u1_k, u2_k, u1_k1, u2_k1,...
     m0, m1, Jg, a, b, k0, k1, mu, c0, c1, T)

 c0 =  1.2*c0;
 k1 = 1.2*k1;
 k0= 1.2*k1;
 
b11_k = cos(x3_k);
b12_k = cos(x3_k);
b21_k = -a;
b22_k = b;


b11_k1 = cos(x3_k1);
b12_k1 = cos(x3_k1);
b21_k1 = -a;
b22_k1 = b;


g1 = x2_k + (3*T/(m0+2*m1))*(-k0*x1_k-k1*x1_k^3-c0*x2_k-c1*(x2_kd^2)*sign(x2_kd)+b11_k*u1_k+b12_k*u2_k) ...
    - (T/(m0+2*m1))*(-k0*x1_k1-k1*x1_k1^3-c0*x2_k1-c1*(x2_kd1^2)*sign(x2_kd1)+b11_k1*u1_k1+b12_k1*u2_k1);

g2 = x4_k + (3*T/(Jg))*(-mu*sign(x4_k)+b21_k*u1_k+b22_k*u2_k)...
    - (T/(Jg))*(-mu*sign(x4_k1)+b21_k1*u1_k1+b22_k1*u2_k1);

y(1) = x1_k + (3*T/2)*x2_k - (T/2)*x2_k1;
y(2) = g1;% + (3*T/2)*b11_k*u1_k + (3*T/2)*b12_k*u2_k;
y(3) = x3_k + (3*T/2)*x4_k - (T/2)*x4_k1;
y(4) = g2;% + (3*T/2)*b21_k*u1_k + (3*T/2)*b22_k*u2_k;
