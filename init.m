%% Denoising ECG signal using ICA, DWT, and PSO. 

%% Initialising
clc;clear all;close all;
fprintf('---- Loading the ECG signal, and computing power line interference ------\n\n');
rng(42);
[signal,Fs]=audioread('Greeting.wav');
z = awgn(signal',5,'measured');              
figure
plot(z);
hold on;
plot(signal);
title('Original signal');

%% PSO
fprintf('--- Setting up PSO parameters ---\n');
nVar=2;            % Number of Decision Variables0
VarSize=[1 nVar];   % Size of Decision Variables Matrix
Var1Min=1;         % Lower Bound of Variables
Var1Max= fix(log2(length(signal)));         % Upper Bound of Variables
%Var1Max=5;
Var2Min=-10;
Var2Max=10;
%PSO Params
MaxIt=70;      % Iterations
nPop=10;        % (Swarm Size)
w=1;            % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=1.1;         % Learning Coefficient
c2=1.6;         % Learning Coefficient
Vel1Max=0.1*(Var1Max-Var1Min);
Vel1Min=-Vel1Max;
Vel2Max=0.1*(Var2Max-Var2Min);
Vel2Min=-Vel2Max;


empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
particle=repmat(empty_particle,nPop,1);
GlobalBest.Cost=-inf;
fprintf('--- Setting initial parameters ---\n');
for i=1:nPop
    particle(i).Position=unifrnd([Var1Min, Var2Min],[Var1Max, Var2Max]);
    particle(i).Position(1)=ceil(particle(i).Position(1));
    particle(i).Velocity=zeros(VarSize);
    particle(i).Cost=CostFunction(signal,z,particle(i).Position);
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    if particle(i).Best.Cost>GlobalBest.Cost
        GlobalBest=particle(i).Best;
    end
    fprintf('--- Completed setting up particle %d ---\n',i);
end
BestCost=zeros(MaxIt,1);

fprintf('--- Starting PSO ---\n');
for it=1:MaxIt    
    for i=1:nPop
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        particle(i).Velocity(1) = max(particle(i).Velocity(1),Vel1Min);
        particle(i).Velocity(1) = min(particle(i).Velocity(1),Vel1Max);
        particle(i).Velocity(2) = max(particle(i).Velocity(2),Vel2Min);
        particle(i).Velocity(2) = min(particle(i).Velocity(2),Vel2Max);
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        if particle(i).Position(1)<Var1Min | particle(i).Position(1)>Var1Max
            particle(i).Velocity(1) = -particle(i).Velocity(1);
        end
        if particle(i).Position(2)<Var2Min | particle(i).Position(2)>Var2Max
            particle(i).Velocity(2) = -particle(i).Velocity(2);
        end
        
        particle(i).Position(1) = ceil(particle(i).Position(1));        
        particle(i).Position(1) = max(particle(i).Position(1),Var1Min);
        particle(i).Position(1) = min(particle(i).Position(1),Var1Max);
        particle(i).Position(2) = max(particle(i).Position(2),Var2Min);
        particle(i).Position(2) = min(particle(i).Position(2),Var2Max);
        particle(i).Cost = CostFunction(signal,z,particle(i).Position);
        
        if particle(i).Cost>particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            if particle(i).Best.Cost>GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end 
        end
        fprintf('--- Iteration %d Particle %d Position %d %d Velocity %d %d Cost %d ---\n',it,i,particle(i).Position(1),particle(i).Position(2),particle(i).Velocity(1), particle(i).Velocity(2),particle(i).Cost);
    end
    BestCost(it)=GlobalBest.Cost;
    fprintf('Iteration  %d Best Cost %d ---\n',it,BestCost(it));
    w=w*wdamp; 
end
BestSol = GlobalBest;
fprintf('--- Completed PSO ---\nComputing denoised signal\n');
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
figure
plot(z);hold on;
plot(mod_sig);
title('Recovered signal');
figure
subplot(1,3,1); specgram(signal,512,Fs); title('True Speech Signal');
subplot(1,3,2); specgram(z,512,Fs); title('Noisy Speech Signal');
subplot(1,3,3); specgram(mod_sig,512,Fs); title('De-noised Speech Signal');

%% Only DWT and iDWT

[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('db13');
[c,ll]=wavedec(z,3,Lo_D,Hi_D);
A=wrcoef('a',c,ll,Lo_R,Hi_R,3);
mod_sig=A;
for i=1:3
    D = wrcoef('d',c,ll,Lo_R,Hi_R,i);
    thr = thselect(D,'minimaxi');
    tD = wthresh(D,'s',thr);
    mod_sig=mod_sig+tD;
end
cost = 20*log10(norm(signal(:)) / norm (signal(:)-mod_sig(:)));
fprintf('\n--- SNR (if only DWT) %d ---\n',cost);

