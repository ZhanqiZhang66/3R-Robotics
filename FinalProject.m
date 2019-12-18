
% q1 = 10;
% q2 = 20;
% q3 = 30;
L1 = 4;
L2 = 3;
L3 = 2;

syms q1 q2 q3 g;
% syms L1 L2 L3 m1 m2 m3



m1 = 20;
m2 = 15;
m3 = 10;


Ic1 = [0 0 0;
    0 0 0;
    0 0 0.5];

Ic2 = [0 0 0;
    0 0 0;
    0 0 0.2];

Ic3 = [0 0 0;
    0 0 0;
    0 0 0.1];


T_01 = [cos(q1) -sin(q1) 0 0;
    sin(q1) cos(q1) 0 0;
    0 0 1 0;
    0 0 0 1];

T_02 = [cos(q1+q2) -sin(q1+q2) 0 L1*cos(q1);
    sin(q1+q2) cos(q1+q2) 0 L1*sin(q1);
    0 0 1 0;
    0 0 0 1];

T_03 = [cos(q1+q2+q3) -sin(q1+q2+q3) 0 L2*cos(q1+q2)+L1*cos(q1);
    sin(q1+q2+q3) cos(q1+q2+q3) 0 L2*sin(q1+q2)+L1*sin(q1);
    0 0 1 0;
    0 0 0 1];

T_0E = [cos(q1+q2+q3) -sin(q1+q2+q3) 0 L3*cos(q1+q2+q3)+L2*cos(q1+q2)+L1*cos(q1);
    sin(q1+q2+q3) cos(q1+q2+q3) 0 L3*sin(q1+q2+q3)+L2*sin(q1+q2)+L1*sin(q1);
    0 0 1 0;
    0 0 0 1];

J0 = [-L3*sin(q1+q2+q3)-L2*sin(q1+q2)-L1*sin(q1) -L3*sin(q1+q2+q3)-L2*sin(q1+q2) -L3*sin(q1+q2+q3);
    L3*cos(q1+q2+q3)+L2*cos(q1+q2)+L1*cos(q1) L3*cos(q1+q2+q3)+L2*cos(q1+q2) L3*cos(q1+q2+q3);
    0 0 0;
    0 0 0;
    0 0 0;
    1 0 0];




% Mass matrix
% Jv1 = [diff(T_01(1:3,1:3) * [L1/2;0;0], q1) diff(T_01(1:3,1:3) * [L1/2;0;0], q2) diff(T_01(1:3,1:3) * [L1/2;0;0], q3)];
% Jv2 = [diff(T_02(1:3,1:3) * [L2/2;0;0], q1) diff(T_02(1:3,1:3) * [L2/2;0;0], q2) diff(T_02(1:3,1:3) * [L2/2;0;0], q3)];
% Jv3 = [diff(T_03(1:3,1:3) * [L3/2;0;0], q1) diff(T_03(1:3,1:3) * [L3/2;0;0], q2) diff(T_01(1:3,1:3) * [L3/2;0;0], q3)];

% Jv1 = [diff(T_01 * [L1/2;0;0;1], q1) diff(T_01 * [L1/2;0;0;1], q2) diff(T_01 * [L1/2;0;0;1], q3)];
% Jv2 = [diff(T_02 * [L2/2;0;0;1], q1) diff(T_02 * [L2/2;0;0;1], q2) diff(T_02 * [L2/2;0;0;1], q3)];
% Jv3 = [diff(T_03 * [L3/2;0;0;1], q1) diff(T_03 * [L3/2;0;0;1], q2) diff(T_01 * [L3/2;0;0;1], q3)];
% Jv1 = Jv1(1:3,1:3);
% Jv2 = Jv2(1:3,1:3);
% Jv3 = Jv3(1:3,1:3);

Pc1 = T_01 * [L1/2;0;0;1];
Pc1 = Pc1(1:3);

Pc2 = T_02 * [L2/2;0;0;1];
Pc2 = Pc2(1:3);

Pc3 = T_03 * [L3/2;0;0;1];
Pc3 = Pc3(1:3);

Jv1 = [diff(Pc1,q1) diff(Pc1,q2) diff(Pc1,q3)];
Jv2 = [diff(Pc2,q1) diff(Pc2,q2) diff(Pc2,q3)];
Jv3 = [diff(Pc3,q1) diff(Pc3,q2) diff(Pc3,q3)];

Jw1 = [0 0 0;
    0 0 0;
    1 0 0];
Jw2 = [0 0 0;
    0 0 0;
    1 1 0];
Jw3 = [0 0 0;
    0 0 0;
    1 1 1];


M = zeros(3);
M = m1*(Jv1.')*Jv1 + m2*(Jv2.')*Jv2 + m3*(Jv3.')*Jv3 + ...
    (Jw1.')*Ic1*Jw1 + (Jw2.')*Ic2*Jw2 + (Jw3.')*Ic3*Jw3;

M = simplify(M);
% V matrix
m111 = diff(M(1,1),q1);
m112 = diff(M(1,1),q2);
m113 = diff(M(1,1),q3);

m121 = diff(M(1,2),q1);
m122 = diff(M(1,2),q2);
m123 = diff(M(1,2),q3);

m131 = diff(M(1,3),q1);
m132 = diff(M(1,3),q2);
m133 = diff(M(1,3),q3);

m211 = diff(M(2,1),q1);
m212 = diff(M(2,1),q2);
m213 = diff(M(2,1),q3);

m221 = diff(M(2,2),q1);
m222 = diff(M(2,2),q2);
m223 = diff(M(2,2),q3);

m231 = diff(M(2,3),q1);
m232 = diff(M(2,3),q2);
m233 = diff(M(2,3),q3);

m311 = diff(M(3,1),q1);
m312 = diff(M(3,1),q2);
m313 = diff(M(3,1),q3);

m321 = diff(M(3,2),q1);
m322 = diff(M(3,2),q2);
m323 = diff(M(3,2),q3);

m331 = diff(M(3,3),q1);
m332 = diff(M(3,3),q2);
m333 = diff(M(3,3),q3);


C = 0.5*[(m111 + m111 - m111), (m122 + m122 - m221), (m133 +m133 - m331);
    (m211 + m211 - m112), (m222 + m222 - m222), (m233 + m233 - m332);
    (m311 + m311 - m113), (m322 + m322 - m223), (m333 + m333 - m333)];

C = simplify(C);

B = [(m112 + m121 - m121), (m113 + m131 - m131), (m123 + m132 - m231);
    (m212 + m221 - m122), (m213 + m231 - m132), (m223 + m232 - m232);
    (m312 + m321 - m123), (m313 + m331 - m133), (m323 + m332 - m233)];
B = simplify(B);

% G matrix
G =  -[(Jv1.') (Jv2.') (Jv3.')] * [m1*[0;-g;0];  m2*[0;-g;0]; m3*[0;-g;0]];

%cartesian
 J01 = [-L3*sin(q1+q2+q3)-L2*sin(q1+q2)-L1*sin(q1) -L3*sin(q1+q2+q3)-L2*sin(q1+q2) -L3*sin(q1+q2+q3);
        L3*cos(q1+q2+q3)+L2*cos(q1+q2)+L1*cos(q1) L3*cos(q1+q2+q3)+L2*cos(q1+q2) L3*cos(q1+q2+q3);
        1 0 0;];



Mx = inv(J01.')* M*inv(J01)



