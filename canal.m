%%%%%%%%%%%%%%%%%%%% Canal %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
             
R = 4;

Tb = 1/R
T=1;
t=0:0.001:T;
f=10;

%generate a signal 
y = sin(2*pi*f*t)
figure
plot(y)

%time shift
Rr = [0 Tb];
v = rand(1,1)*range(Rr)+min(Rr)
A = zeros(1,round(v));
C = horzcat(A,y)
figure
plot(C);

%amortissement 
Rt = [0 0.99];
alphaN = rand(1,1)*range(Rt)+min(Rt)
D = (C.* alphaN);
figure
plot(D)

%somme les signaux de chaque canal
E = vertcat(D,D,D);
S = sum(E);
figure 
plot(S)

%calcul de l'énergie de l'AWGN dans la BP
No = Eb/EbNo;
%Densite spectrale de puissance [Watts/Hz]
E_AWGN = N0/2/T_b*2*N; 

%ajouter bruit blanc gausien
snr = 10; %Eb/No
%G = awgn(S,snr)
%figure
%plot(G)