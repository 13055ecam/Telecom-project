%%%%%%%%%%%Rajouter des bits de synchronisation comme 101010%%%%%%%%%%%%%%
%%%%%On code en PAM binaire car entre -1 et 1 c est plus facile a
%%%%%distinguer qu entre 0 et 1

%Bit rate in bits/s
R = 10;
%Seconds per bit = symbol duration
Tb = 1/R;
% M is the number of bits = number of symbols in each weft (trame)
M = 10;
k = 0:M-1;
%Message xn(k)
xn = 2*randi([0 1], 1, M) - 1;

an = xn;

%Number of material resources available
N = 2;
n = 0:N-1;

%Amplitude
A = 10;

beta = (4*N)-2;
Tn = Tb/beta;

wn = 2*pi*2*n/Tb;

rf = 0;
span = 10;

h2 = rcosdesign(rf,span,beta,'sqrt');
% fvtool(h2,'impulse')

codes_fir = upfirdn(an, h2, beta);
sizeupfirdn = size(codes_fir, 2);

porteuse = cos(wn'.*linspace(0, Tn*sizeupfirdn, sizeupfirdn));

s = porteuse.*codes_fir;

% stem(s(1,:))
sfft = fft(s');
length(linspace(0, sizeupfirdn*Tn, sizeupfirdn))
length(s)
figure
plot((ones(1,2)'*linspace(0, sizeupfirdn*Tn, sizeupfirdn))', s')
subplot(2,1,1)
stem((ones(1,2)'*linspace(0, sizeupfirdn*Tn, sizeupfirdn))', s')
grid

% figure

% plot((ones(1,2)'*linspace(0, (1/Tn)-1, sizeupfirdn))', 20*log10(abs(sfft)'/(M*beta)))
subplot(2,1,2)
stem((ones(1,2)'*linspace(0, (1/Tn)-1, sizeupfirdn))', abs(sfft)/(M*beta))
grid

% sHighb = ifft(sHighFFTb);
% stem(linspace(0, sizeupfirdn(2)*Tn, sizeupfirdn(2)), sHighb)
% grid
% subplot(2,2,4)
% stem(linspace(0, 1/Tn, sizeupfirdn(2)), abs(sHighFFTb)/(M*beta))
% grid

%Calculer le bon facteur d amplitude permettant que le signal de tension
%corresponde a la puissance qu on veut (quelques mW), on utilise pour ca
%l impedance caracteristique d un cable. Prendre 75ohms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% porteuse = cos(linspace(-L*Tb, L*Tb, sizeupfirdn(2)).*wn');
% porteuse = cos(fn'/(Tn*sizeupfirdn(2)).*linspace(0, 2*pi*Tn*sizeupfirdn(2), sizeupfirdn(2)));
% L = 10;
% temps = ((-L*Tb):Tn:(L*Tb));
% 
% upsampled = upsample(an, beta);
% 
% rf = 0.25;
% span = 10;
% sps = 2;
% h2 = rcosdesign(rf,span,sps,'sqrt');
% 
% s = filter(h2, 1, upsampled);
% 
% % pn = h2.*cos(wn'*temps);
% pn_fft = fft(pn);
% 
% % sn = conv(upsampled, pn);
% 
% 
% %filtre en cos sureleve
% % h = rcosdesign(0.25,6,4, 'sqrt');
% % mx = max(abs(h-rcosdesign(0.25,6,4,'sqrt')))
% % fvtool(h,'Analysis','impulse')
% 
% rf = 0;
% span = 10;
% sps = 20;
% 
% h1 = rcosdesign(rf,span,sps,'normal');
% fvtool(h1,'impulse');

% h2 = rcosdesign(rf,span,sps,'sqrt');
% fvtool(h2,'impulse');
% 
% h3 = conv(h2,h2);
% p2 = ceil(length(h3)/2);
% m2 = ceil(p2-length(h1)/2);
% M2 = floor(p2+length(h1)/2);
% ct = h3(m2:M2);
% 
% stem([h1/max(abs(h1));ct/max(abs(ct))]','filled')
% xlabel('Samples')
% ylabel('Normalized amplitude')
% legend('h1','h2 * h2')
% 
% b = rcosdesign(rf, span, sps);
% d = 2*randi([0 1], 100, 1) - 1;
% x = upfirdn(d, b, sps);
% r = x + randn(size(x))*0.01;
% y = upfirdn(r, b, 1, sps);
% figure
% % stem(y)
% t = randi([0 1], 1, length(wn))
% wn'
% t*wn'
% length(wn)
% % t = -2*pi:0.01:2*pi;
% cos(wn'*t)
% pn = h2.*cos(wn'*temps);
% figure
% stem(pn)
