
clc; clear all; close all;

addpath(genpath('C:\Users\9301446\Downloads\YALMIP-master'))
warning off;

% INITIAL SYSTEM PARAMETERS
T = 5e-3;
tf = .5;
stp = ceil(tf/T);
t = 0: T: (stp-1)*T;

m0 = 1;
m1 = 0.1;
k0 = 0.2;
k1 = 0.1;
c0 = 0.3;
c1 = 0.05;
mu = 0.3;

l = 1;
a = l/3;
b = l/4;

Jg = (1/12)*m0*(4*l^2)+2*m1*l^2;


for i = 1:stp
    u(i) = 0;
end

x = zeros(3,4);
u1 = zeros(3,1);
u2 = zeros(3,1);

ax1 = 1.5;
wx1 = 2*pi;
ax3 = .5;
wx3 = 1.5*pi;


stpf = 10;
for i = 1:stp+2
    i*T
    x1d(i) = ax1*cos(wx1*i*T);
    x3d(i) = ax3*cos(wx3*i*T);
end

lambda = 0.1;
alpha1 = 0.9;
alpha2 = 0.9;


umax1 = 100;
umax2 = 100;

rho = .01;

k=3

md = 5;
Md = 10;
k=3;

x1max = 2;
x2max = 1;

R = 1e-5;
%%
for k = 3:stp
    %%    SAMPLING PLANT
    t_k = k*T
    
    if k>md+1
        x2_kd = x(k-md,2);
        x2_kd1 = x(k-md-1,2);
    else
        x2_kd = 0;
        x2_kd1 = 0;
    end
    
    x(k,:) = syseq (x(k-1,1), x(k-1,2), x(k-1,3), x(k-1,4),  ...
        x(k-2,1), x(k-2,2), x(k-2,3), x(k-2,4),...
        x2_kd, x2_kd1,...
        u1(k-1), u2(k-1), u1(k-2), u2(k-2),...
        m0, m1, Jg, a,b, k0, k1, mu, c0, c1, T);
    
    
    x1_k = x(k,1);
    x2_k = x(k,2);
    x3_k = x(k,3);
    x4_k = x(k,4);
    
    x1_k1 = x(k-1,1);
    x2_k1 = x(k-1,2);
    x3_k1 = x(k-1,3);
    x4_k1 = x(k-1,4);
    
    u1_k1 = u1(k-1);
    u2_k1 = u2(k-1);
    
    b11_k = cos(x3_k);
    b12_k = cos(x3_k);
    b21_k = -a;
    b22_k = b;
    
    
    b11_k1 = cos(x3_k1);
    b12_k1 = cos(x3_k1);
    b21_k1 = -a;
    b22_k1 = b;
    
    x1_f1 = x1_k + (3*T/2)*x2_k - (T/2)*x2_k1;
    x3_f1 = x3_k + (3*T/2)*x4_k - (T/2)*x4_k1;
    
  for  i = 1:Md
    if k>Md+1
        x2_kd = x(k-i+1,2);
        x2_kd1 = x(k-i,2);
        
    else
        x2_kd = 0;
        x2_kd1 = 0;
    end
    
    g1d0(i) = (3*T/(m0+2*m1))*(-c1*(x2_kd^2)*sign(x2_kd)) ...
        - (T/(m0+2*m1))*(-c1*(x2_kd1^2)*sign(x2_kd1));
  end
    
      g1d = max(g1d0);
    
      
    g10 = x2_k + (3*T/(m0+2*m1))*(-k0*x1_k-k1*x1_k^3-c0*x2_k) ...
        - (T/(m0+2*m1))*(-k0*x1_k1-k1*x1_k1^3-c0*x2_k1+b11_k1*u1_k1+b12_k1*u2_k1);
    
    g2 = x4_k + (3*T/(Jg))*(-mu*sign(x4_k))...
        - (T/(Jg))*(-mu*sign(x4_k1)+b21_k1*u1_k1+b22_k1*u2_k1);
    
     g1v = g10 + g1d0;
    
    g1 = g1v(1);
    for i=1:max(size(g1v))
        if abs(g1v(i))>abs(g1)
    g1 = g1v(i);
        end
    end
    
    s1(k) = x1_f1 - x1d(k+1) + lambda*(x1_k - x1d(k));
    kesi1 = x1_f1 + (3*T/2)*g1 - (T/2)*x2_k - x1d(k+2) + lambda*x1_f1 - lambda*x1d(k+1);
    
    s2(k) = x3_f1 - x3d(k+1) + lambda*(x3_k - x3d(k));
    kesi2 = x3_f1 + (3*T/2)*g2 - (T/2)*x4_k - x3d(k+2) + lambda*x3_f1 - lambda*x3d(k+1);
    
    av = [kesi1 ; kesi2];
    %%
    b011 = (9*T^2/4)*b11_k;
    b012 = (9*T^2/4)*b12_k;
    b021 = (9*T^2/4)*b21_k;
    b022 = (9*T^2/4)*b22_k;
    
    Bmat = [b011 b012 ; b021 b022];
    
    Bmatinv  = inv(Bmat);
    
    Bconst = [ 1 , 0 ;
        0 , 1;];
    %         (3*T/2)*b11 , (3*T/2)*b12 ;
    %         (3*T/2)*b21 , (3*T/2)*b22 ] ;
    
    qconsth = [ umax1 ; umax2 ;];% accmax1 - h1 + h3 ; accmax1 - h2 + h4];
    qconstl = - [ umax1 ; umax2 ;];% -accmax1 - h1 + h3 ; -accmax1 - h2 + h4];
    
    u10 = sdpvar(1);
    u20 = sdpvar(1);
    

    
    
    %Sv0 = [s10 ; s20];
    Sv = [s1(k) ; s2(k)];
    
    P = sdpvar(2);
    M = sdpvar(2);
    X0 = sdpvar(2);
    Sv0 = M*Sv;
    
        s10 = Sv0(1);
    s20 = Sv0(2);
    numax = [0; 0];
    for i = 0:umax1
        for j = 0:umax2
            umax0 = abs(Bmat*[i ; j]);
            
            if umax0(1)>numax(1)
                numax(1) = umax0(1);
            end
            if umax0(2)>numax(2)
                numax(2) = umax0(2);
            end
            
        end
    end
    
c2 = 1;
    
      for i = 0:.1:x1max
            sfeas10(c2) = abs(i - x1d(k) + lambda*(x(k-1,1)-x1d(k-1)));
            c2= c2+1;
      end
      c2 = 1;
      
              for j = 0:.1:x2max

            sfeas20(c2) = abs(j - x3d(k) + lambda*(x(k-1,3)-x3d(k-1)));
                        c2= c2+1;
              end
        
      sfeas1 = max(sfeas10);
      sfeas2 = max(sfeas20);
    
      sfeas1 = min(sfeas1,abs(s1(k)));
      sfeas2 = min(sfeas2,abs(s2(k)));
     
  Umax = [umax1; umax2];
%     
%     numax = abs(Bmat*Umax); 
%     
    %     if (abs(s1(k))>rho) && (abs(s2(k))>rho)
    
    %         F = [(b011*u10 + b012*u20)<0.99*abs(s1(k))-kesi1sup, (b011*u10 + b012*u20)>-0.99*abs(s1(k))-kesi1inf , ...
    %             (b021*u10 + b022*u20)<0.99*abs(s2(k))-kesi2, (b021*u10 + b022*u20)>-0.99*abs(s2(k))-kesi2 , ...
    %             Bconst*[u10 ; u20]<qconsth , Bconst*[u10 ; u20]>qconstl];
    %
    F = [Sv0.'*(P+Bmatinv.'*R*Bmatinv)*Sv0-2*Sv0.'*Bmatinv.'*R*Bmatinv*av-Sv.'*Sv+av.'*Bmatinv.'*R*Bmatinv*av<0 , abs(s10)<sfeas1, abs(s20)<sfeas2,...
        Sv0.'*P*Sv0<1,Sv.'*P*Sv<1,P(1,1)>0,P(1,1)*P(2,2)-P(1,2)^2>0];%...
%         ,abs(Bmat\(Sv0-av))<Umax];
    f = [Sv0.'*Sv0+(Sv0.'*(Bmatinv.'*Bmatinv)*Sv0-2*Sv0.'*Bmatinv.'*Bmatinv*av+av.'*Bmatinv.'*Bmatinv*av)*R]
    %     else
    %
    %         F = [(b011*u10 + b012*u20)<0.99*rho-kesi1, (b011*u10 + b012*u20)>-0.99*rho-kesi1 , ...
    %             (b021*u10 + b022*u20)<0.99*rho-kesi2, (b021*u10 + b022*u20)>-0.99*rho-kesi2 , ...
    %             Bconst*[u10 ; u20]<qconsth , Bconst*[u10 ; u20]>qconstl ];
    %
    %     end
    
    %  f = ( b11*u10 + b12*u20 + h1(k) - h3(k))^2+( b21*u10 + b22*u20 + h2(k) - h4(k))^2;
    
    sol=solvesdp(F,f);
    
    %      if (abs(double(u10))<umax1) && (abs(double(u20))<umax2)
    %         u1(k) = double(u10);
    %         u2(k) = double(u20);
    %
    mode(k) = 1;
    %     else
    %
    %
    %         %%
    B = [b011 b012 ; b021 b022];
    %
    %     qh = [abs(s1(k))-kesi1 ; abs(s2(k))-kesi2];
    %     ql = [-abs(s1(k))-kesi1 ; -abs(s2(k))-kesi2];
    %
    %     alpha = [alpha1 0 ; 0 alpha2];
    %     qm = (eye(2)-alpha)*ql + alpha*qh;
    %
    %
    %     uv = B\qm;
    %
    Sv0 = double(Sv0);
    %Sv0 = [0 ; 0];
    uvp = B\(Sv0-av);
    %%
    u1(k) = double(uvp(1));
    u2(k) = double(uvp(2));
    %
        thres1 = umax1;
        thres2 = umax2;
    
        if (abs(u1(k)))>thres1
            u1(k) = thres1*sign(u1(k));
        end
    
        if (abs(u2(k)))>thres2
            u2(k) = thres2*sign(u2(k));
        end
    
    mode(k) = 2;
    
    
    
    
end
%%
close all;

figure; plot(t,x(:,1)); hold all; plot(t,x1d(1:max(size(t)))); legend('System Response' , 'Reference Signal');  xlabel('Time [s]'); ylabel('X_1 [m]');
figure; plot(t,x(:,2)); xlabel('Time [s]'); ylabel('X_2 [m/s]');

figure; plot(t,x(:,3)); hold all; plot(t,x3d(1:max(size(t)))); legend('System Response' , 'Reference Signal');  xlabel('Time [s]'); ylabel('X_3 [m]');
figure; plot(t,x(:,4)); xlabel('Time [s]'); ylabel('X_4 [m/s]');

figure; plot(t(2:end),u1(2:end)); xlabel('Time [s]'); ylabel('u_1 [m]');
figure; plot(t(2:end),u2(2:end)); xlabel('Time [s]'); ylabel('u_2 [m]');
