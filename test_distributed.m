clear all; clc;

J1 = rand(5,3);
J2 = rand(7,6);
D1 = rand(5,2);
D2 = rand(7,2);


J = [J1, zeros(size(J1,1), size(J2,2)), D1; zeros(size(J2,1), size(J1,2)) J2 D2];

N = J'*J;

R1 = D1'*D1 - D1'*J1*inv(J1'*J1)*J1'*D1
R2 = D2'*D2 - D2'*J2*inv(J2'*J2)*J2'*D2

