%% Denoising ECG signal DWT, and PSO. 

%% Initialising
clc;clear all;close all;
fprintf('---- Loading the signal, and computing signal with power line interference ------\n\n');
rng(42);
% [signal,Fs]=audioread('Greeting.wav');
Fs=360;
load('105m.mat');
signal = val(1,:);
signal=(signal-0)/200;
%z = awgn(signal',5,'measured');              
x=0.1*(max(signal)-min(signal));
t=(0:length(signal)-1)/Fs;
z=signal+0.2*sin(2*pi*50*t);
figure
subplot(2,1,1);
plot(z);
title('Signal with Power line interference');
subplot(2,1,2);
plot(signal);
title('Signal');

%% PSO
BestSol = PSO(signal, z, 100, 100, 1);
%% Calculating final signal
fprintf('Level %d Threshold %d\n',BestSol.Position(1), BestSol.Position(2));
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('db13');
[c,ll]=wavedec(z,BestSol.Position(1),Lo_D,Hi_D);
A=wrcoef('a',c,ll,Lo_R,Hi_R,BestSol.Position(1));
mod_sig=A;
for i=1:BestSol.Position(1)
    D = wrcoef('d',c,ll,Lo_R,Hi_R,i);
    tD = wthresh(D,'s',BestSol.Position(2));
    mod_sig=mod_sig+tD;
end
cost = 20*log10(norm(signal(:)) / norm (signal(:)-mod_sig(:)));
fprintf('--- SNR %d ---\n',cost);
%% Calculating signal using inbuilt method
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('db13');
[c,ll]=wavedec(z,3,Lo_D,Hi_D);
A=wrcoef('a',c,ll,Lo_R,Hi_R,3);
mod_sig2=A;
for i=1:3
    D = wrcoef('d',c,ll,Lo_R,Hi_R,i);
    thr = thselect(D,'minimaxi');
    tD = wthresh(D,'s',thr);
    mod_sig2=mod_sig2+tD;
end
cost2 = 20*log10(norm(signal(:)) / norm (signal(:)-mod_sig2(:)));
fprintf('\n--- SNR (if only DWT) %d ---\n',cost2);
%% Plotting
figure
subplot(2,1,1);
plot(mod_sig2);
title('Using ThSelect and Level3');
subplot(2,1,2);
plot(mod_sig);
title('Using PSO');
% figure
% subplot(1,4,1); specgram(signal,512,Fs); title('True Speech Signal');
% subplot(1,4,2); specgram(z,512,Fs); title('Noisy Speech Signal');
% subplot(1,4,3); specgram(mod_sig,512,Fs); title('De-noised Speech Signal');
% subplot(1,4,4); specgram(mod_sig2,512,Fs); title('De-noised Speech Signal Using Inbuilt');