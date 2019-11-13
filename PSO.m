function position = PSO(signal, z, popsize, iterations, printflag)
%% Function to calculate the position required for maximum SNR
%%
if printflag == 1
    fprintf('--- Setting up PSO parameters ---\n');
end
nVar=2;            % Number of Decision Variables0
VarSize=[1 nVar];   % Size of Decision Variables Matrix
Var1Min=1;         % Lower Bound of Variables
Var1Max= fix(log2(length(signal)));         % Upper Bound of Variables
%Var1Max=5;
Var2Min=-(max(z)*2);
Var2Max=-Var2Min;
%PSO Params
MaxIt=iterations;      % Iterations
nPop=popsize;        % (Swarm Size)
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
if printflag == 1
    fprintf('--- Setting initial parameters ---\n');
end
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
    if printflag == 1
        fprintf('--- Completed setting up particle %d ---\n',i);
    end
end
BestCost=zeros(MaxIt,1);

if printflag == 1
    fprintf('--- Starting PSO ---\n');
end
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
        
        if particle(i).Position(1)<Var1Min || particle(i).Position(1)>Var1Max
            particle(i).Velocity(1) = -particle(i).Velocity(1);
        end
        if particle(i).Position(2)<Var2Min || particle(i).Position(2)>Var2Max
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
        if printflag == 1
            fprintf('--- Iteration %d Particle %d Position %d %d Velocity %d %d Cost %d ---\n',it,i,particle(i).Position(1),particle(i).Position(2),particle(i).Velocity(1), particle(i).Velocity(2),particle(i).Cost);
        end
    end
    BestCost(it)=GlobalBest.Cost;
    if printflag == 1
        fprintf('Iteration  %d Best Cost %d ---\n',it,BestCost(it));
    end
    w=w*wdamp; 
end
position = GlobalBest;
end