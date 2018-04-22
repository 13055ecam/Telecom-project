%%%%%%%%%%%%%%%%%%%% Canal %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

%parameters of a signal 
R = 4;
Tb = 1/R;
T=1;
t=0:0.001:T;
f=10;

%generate a signal 
y = sin(2*pi*f*t);
figure
plot(y)

%time shift
Rr = [0 Tb]; %delais plus petit que Tb
v = rand(1,1)*range(Rr)+min(Rr); %random delais
Add_zero = zeros(1,round(v)); %on ajoute "v" zéros pour décaller le signal 
shifted_data = horzcat(Add_zero,y);%les zeros devant puis le signal après 
figure
plot(shifted_data);

%amortissement 
Rt = [0 0.99]; %facteur d’affaiblissement plus petit que 1
alphaN = rand(1,1)*range(Rt)+min(Rt) %random facteur d'affaiblissement
damped_data = (shifted_data.* alphaN);
figure
plot(damped_data)

%somme les signaux de chaque canal
E = vertcat(D,D,D);
S = sum(E);
figure 
plot(S)

%ajout du bruit blanc gausien n(t) No/2
snr = 10; %Eb/No
G = awgn(S,snr)
figure
plot(G)
