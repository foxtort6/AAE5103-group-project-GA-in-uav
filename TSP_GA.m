
%UAV mission points shortest paths 

clear all;clc
x = randi([-100 100], 1, 30);% Randomly set 30 integers between -100 to 100 to represent the mission points coordinations
y = randi([-100 100], 1, 30);
E = (1:30);


missions=[E' x' y']; %use those numbers to denote mission point possition

scatter(missions(:,2),missions(:,3)) %draw mission points on the figure

%initial parameters
PC=0.8; %crossover probability
PM=0.2;  %mutation probability
iter=2500;   %iteration
pop=200; %population size
n=30;   %chromosome length

initial_pop= zeros(pop,n);
for i =1:pop
    initial_pop(i,:)=randperm(n); %initial population
end

%major program of GA
for genetations=1:iter

%Encoding and Ranking
paths=zeros(pop,n+1);   %Record 31 paths(including last one) of iterations
distance=zeros(pop,1);  %Record every distances of iterations
fitness=zeros(pop,1);   %Record every fitness of iterations

for i=1:pop   
    paths(i,:)=[initial_pop(i,:),initial_pop(i,1)]; %go back to origin
    
    for j =1:n  %Calculating the Euclidean Distance between every mission points
        distance(i)=distance(i)+sqrt((missions(paths(i,j),2)-missions(paths(i,j+1),2))^2+(missions(paths(i,j),3)-missions(paths(i,j+1),3))^2);
    end
    
    fitness(i)=1/distance(i); %fitness calculating
end

[a,b]=max(fitness);  % a record maximum of fitness b record the row
[c,d]=min(distance);  % c record minimum of the fitness d record the row 

%Selection
selection_pop=zeros(pop,n);

FIT1=(fitness-min(fitness))/(max(fitness)-min(fitness)); %transferring all the fitness between 0 to 1
FIT2=FIT1/sum(FIT1);
FIT3=cumsum(FIT2);

BOX=sort(rand(pop,1));  %Roulette strategy: random generate 200 numbers from 0 to 1 in a box
i=1;j=1;

while i<=pop                                                               
    if BOX(i)<FIT3(j)
        selection_pop(i,:)=initial_pop(j,:);   %Compare the fit with the random number in Box
        i=i+1;
    else
        j=j+1;
    end
end

%Crossover
crossover_pop=selection_pop;

for i=1:2:pop-1  %Crossover between two adjacent chromosomes.
    if rand<PC
        store=crossover_pop(i,10:15);  %to store the crossover position 
        
        for k=1:length(store)
            f(k)=find(crossover_pop(i+1,:)==store(k)); %find the correlated position in matching chromosome
        end
        
        f=sort(f);  %Sort to record the sequence of the position of genes.
        
        for p=1:length(store)
            crossover_pop(i,p+9)=crossover_pop(i+1,f(p));  
            crossover_pop(i+1,f(p))=store(p);
        end
        
    end
end

%Mutation
mutation_pop=crossover_pop;

for i=1:pop
    if rand<PM
        r=randperm(30);
        r1=r(1);
        r2=r(2);
        row=mutation_pop(i,:);
        store1=find(row==r1);
        store2=find(row==r2);
        mutation_pop(i,store1)=r2;
        mutation_pop(i,store2)=r1;
    end
end
    
%Next generation
mutation_pop(end,:)=initial_pop(b,:); %find best iteration and store them

shortest_distance(genetations)=distance(b);
best_path=paths(b,:);
initial_pop=mutation_pop;
 
drawnow
scatter(missions(:,2),missions(:,3))
hold on
plot(missions(best_path(1:n+1),2),missions(best_path(1:n+1),3),'b')  %figure of plot
title(['No ' num2str(genetations) ' iteration']);
hold off

end

shortest_distance(genetations)  %shortest distance
best_path        %best route

plot(1:iter,shortest_distance);  %draw iteration plot
title 'Distances - Iterations'
figure
plot(missions(best_path(1:n+1),2),missions(best_path(1:n+1),3),'b')  %distance map