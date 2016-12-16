A= [ 
1, 2.593456468536017, 10, 5;
1/2, 2.660722422908167, 20, 10;
1/4, 2.71063630399989, 40, 20;
1/8, 2.699798002375426, 80, 40;
1/16, 2.694352105835493, 160, 80;
1/32, 2.726326573415554, 320, 160;
1/64, 2.718113796848138, 640, 320;
];
B= [
1, 2.687682606413367, 32, 0;
0.1, 2.687682606413367, 32, 0;
0.01, 2.70407465798737, 40, 0;
0.001, 2.716904615962541, 104, 0;
0.0001, 2.718139681269396, 304, 0;
1e-05, 2.718267502076198, 975, 0;
];


X=A(:,2);
Y=A(:,3);
Z=B(:,2);
T=B(:,3);
X=abs(X-exp(1))
Z=abs(Z-exp(1))
plot(X,Y,Z,T);
legend('fixed timesteps','adapted timesteps');