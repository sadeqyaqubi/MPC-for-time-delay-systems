clc; clear all; close all;
%% PMIN CALCULATION
umax = [2.7 ; 2.6 ];
umin = [2.3 ; 2.4 ];
u = [2.4 ; 2.5 ];

mu = min((umax-u),(u - umin));
c1 = 1;
c2 = 1;
alpha = 0.2;
Psi = diag([alpha alpha]);


pf = .5;
ps = 0.05;
xf = 1;
xs = 0.05;

for p11= 0:ps:pf
    p11;
    for p12 = 0:ps:pf;
        for p22 = 0:ps:pf;
            P = [p11 p12; p12 p22];
            checfac = P - Psi.'*P*Psi;
            if (checfac>0)
                p11can(c1) = p11;
                p12can(c1) = p12;
                p22can(c1) = p22;
                c1 = c1+1
            end
        end
    end
end

for i = 1:c1-1
    i
    for x11= 0:ps:mu(1)
        for x12 = 0:ps:pf;
            for x22 = 0:ps:mu(2);
                X = [x11 x12; x12 x22];
                Pmat = [p11can(i) p12can(i); p12can(i) p22can(i)];
                Xmat = [X , Psi; Psi.' , Pmat];
                eigvec = eig(Xmat);
                ceig(c2) = 1;
                
                p11can0(c2) = p11can(i);
                p12can0(c2) = p12can(i);
                p22can0(c2) = p22can(i);
                
                for j = 1:max(size(eigvec))
                    if (real(eigvec(j)))<0
                        ceig(c2) = 0;
                    end
                end
                
                c2 = c2+1;
            end
        end
    end
end

c3 = 1;
for i = 1:c2-1
    if (ceig(i)==1)
        p11can1(c3) = p11can0(i);
        p12can1(c3) = p12can0(i);
        p22can1(c3) = p22can0(i);
        c3 = c3+1;
    end
end
c3
if (c3>1)
    Pmin = [p11can1(1) p12can1(1); p12can1(1) p22can1(1)];
    for i = 1:c3-1
        Pmincan = [p11can1(i) p12can1(i); p12can1(i) p22can1(i)];
        if (trace(Pmincan)<trace(Pmin))
            Pmin = [p11can1(i) p12can1(i); p12can1(i) p22can1(i)];
        end
    end
    
    Pmin
end
A = P0;

%%
A = double(P)
theta = 0:0.1:2.1*pi;
n = size(theta,2); 
x = zeros(theta,1);
y = zeros(theta,1);
for i = 1:n,
angle = theta(i);
c = cos(angle);
s = sin(angle);
X = [ c, s ];
r = sqrt(1/abs(X*A*X.')); 
x(i) = r*c;
y(i) = r*s;
end;
figure; plot(x,y); xlabel('S_1'); ylabel('S_2'); 