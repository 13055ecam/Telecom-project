fc = 300;
fs = 1000;
filter_order = 6;
cutoff_freq = fc/(fs/2)
[b,a] = butter(filter_order,cutoff_freq);
freqz(b,a)

dataIn = randn(1000,1);
dataOut = filter(b,a,dataIn);

[A,B,C,D] = butter(10,[500 560]/750);
d = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',500,'HalfPowerFrequency2',560, ...
    'SampleRate',1500);

sos = ss2sos(A,B,C,D);
fvt = fvtool(sos,d,'Fs',1500);
legend(fvt,'butter','designfilt')



received_signal = randn(1000,1);

dataOut = filter(b,a,received_signal);

