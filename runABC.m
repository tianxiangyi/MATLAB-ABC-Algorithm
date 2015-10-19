%%%%%ARTIFICIAL BEE COLONY ALGORITHM%%%%

%Artificial Bee Colony Algorithm was developed by Dervis Karaboga in 2005 
%by simulating the foraging behaviour of bees.

%Copyright ? 2008 Erciyes University, Intelligent Systems Research Group, The Dept. of Computer Engineering

%Contact:
%Dervis Karaboga (karaboga@erciyes.edu.tr )
%Bahriye Basturk Akay (bahriye@erciyes.edu.tr)


clear all
close all
clc



% Set ABC Control Parameters
ABCOpts = struct( 'ColonySize',  100, ...   % Number of Employed Bees+ Number of Onlooker Bees 
    'MaxCycles', 1000,...   % Maximum cycle number in order to terminate the algorithm
    'ErrGoal',   1e-20, ...  % Error goal in order to terminate the algorithm (not used in the code in current version)
    'Dim',       5, ... % Number of parameters of the objective function   
    'Limit',   50, ... % Control paramter in order to abandone the food source 
    'lb',  -600, ... % Lower bound of the parameters to be optimized
    'ub',  600, ... %Upper bound of the parameters to be optimized
    'ObjFun' , 'griewank', ... %Write the name of the objective function you want to minimize
    'RunTime',1); % Number of the runs 



GlobalMins=zeros(ABCOpts.RunTime,ABCOpts.MaxCycles);

for r=1:ABCOpts.RunTime
    
% Initialise population
Range = repmat((ABCOpts.ub-ABCOpts.lb),[ABCOpts.ColonySize ABCOpts.Dim]);
Lower = repmat(ABCOpts.lb, [ABCOpts.ColonySize ABCOpts.Dim]);
Colony = rand(ABCOpts.ColonySize,ABCOpts.Dim) .* Range + Lower;

Employed=Colony(1:(ABCOpts.ColonySize/2),:);


%evaluate and calculate fitness
ObjEmp=feval(ABCOpts.ObjFun,Employed);
FitEmp=calculateFitness(ObjEmp);

%set initial values of Bas
Bas=zeros(1,(ABCOpts.ColonySize/2));


GlobalMin=ObjEmp(find(ObjEmp==min(ObjEmp),end));
GlobalParams=Employed(find(ObjEmp==min(ObjEmp),end),:);

Cycle=1;
while ((Cycle <= ABCOpts.MaxCycles)),
    
    %%%%% Employed phase
    Employed2=Employed;
    for i=1:ABCOpts.ColonySize/2
        Param2Change=fix(rand*ABCOpts.Dim)+1;
        neighbour=fix(rand*(ABCOpts.ColonySize/2))+1;
            while(neighbour==i)
                neighbour=fix(rand*(ABCOpts.ColonySize/2))+1;
            end;
        Employed2(i,Param2Change)=Employed(i,Param2Change)+(Employed(i,Param2Change)-Employed(neighbour,Param2Change))*(rand-0.5)*2;
         if (Employed2(i,Param2Change)<ABCOpts.lb)
             Employed2(i,Param2Change)=ABCOpts.lb;
         end;
        if (Employed2(i,Param2Change)>ABCOpts.ub)
            Employed2(i,Param2Change)=ABCOpts.ub;
        end;
        
    end;  %for end

    ObjEmp2=feval(ABCOpts.ObjFun,Employed2);
    FitEmp2=calculateFitness(ObjEmp2);
    [Employed ObjEmp FitEmp Bas]=GreedySelection(Employed,Employed2,ObjEmp,ObjEmp2,FitEmp,FitEmp2,Bas,ABCOpts);
    
    %Normalize
    NormFit=FitEmp/sum(FitEmp);
    
    %%% Onlooker phase  
Employed2=Employed;
i=1;
t=0;
while(t<ABCOpts.ColonySize/2)
    if(rand<NormFit(i))
        t=t+1;
        Param2Change=fix(rand*ABCOpts.Dim)+1;
        neighbour=fix(rand*(ABCOpts.ColonySize/2))+1;
            while(neighbour==i)
                neighbour=fix(rand*(ABCOpts.ColonySize/2))+1;
            end;
         Employed2(i,:)=Employed(i,:);
         Employed2(i,Param2Change)=Employed(i,Param2Change)+(Employed(i,Param2Change)-Employed(neighbour,Param2Change))*(rand-0.5)*2;%find a new value in the neighborhood
         if (Employed2(i,Param2Change)<ABCOpts.lb)
             Employed2(i,Param2Change)=ABCOpts.lb;
         end;
        if (Employed2(i,Param2Change)>ABCOpts.ub)
            Employed2(i,Param2Change)=ABCOpts.ub;
         end;
    ObjEmp2=feval(ABCOpts.ObjFun,Employed2);
    FitEmp2=calculateFitness(ObjEmp2);
    [Employed ObjEmp FitEmp Bas]=GreedySelection(Employed,Employed2,ObjEmp,ObjEmp2,FitEmp,FitEmp2,Bas,ABCOpts,i);
   
   end; % end if
    
   i=i+1;
   if (i==(ABCOpts.ColonySize/2)+1) 
       i=1;
   end;   
end; %end while
    
    
    %%%Memorize Best
 CycleBestIndex=find(FitEmp==max(FitEmp));
 CycleBestIndex=CycleBestIndex(end);
 CycleBestParams=Employed(CycleBestIndex,:);
 CycleMin=ObjEmp(CycleBestIndex);
 
 if CycleMin<GlobalMin 
       GlobalMin=CycleMin;
       GlobalParams=CycleBestParams;
 end
 
 GlobalMins(r,Cycle)=GlobalMin;
 
 %% Scout phase
 ind=find(Bas==max(Bas));
ind=ind(end);
if (Bas(ind)>ABCOpts.Limit)
Bas(ind)=0;
Employed(ind,:)=(ABCOpts.ub-ABCOpts.lb)*(0.5-rand(1,ABCOpts.Dim));%*2+ABCOpts.lb;
%message=strcat('burada',num2str(ind))
end;
ObjEmp=feval(ABCOpts.ObjFun,Employed);
FitEmp=calculateFitness(ObjEmp);
    


  %  fprintf('Cycle=%d ObjVal=%g\n',Cycle,GlobalMin);
    
    Cycle=Cycle+1;

end % End of ABC

end; %end of runs
if ABCOpts.RunTime==1
    semilogy(GlobalMins);
else
    semilogy(mean(GlobalMins),'b');
end
title('Mean of Best function values');
xlabel('cycles');
ylabel('error');
fprintf('Mean =%g Std=%g\n',mean(GlobalMins(:,end)),std(GlobalMins(:,end)));
