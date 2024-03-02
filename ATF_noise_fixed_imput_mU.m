clear; clc;
rng(1)

%% Simulacion sistema ATF  
total_time = 50000;

t = 0;      
time_points = zeros(1, 1000000);  

% Initial conditions
W = 0;
Y = 0;
U = 0;
C = 0;
X = [W, Y, U, C]; 

X_results = zeros(100000, 4);


% Kinetic parameters
g = 0.0004;
gU = 0.0004;
gW = 0.0004;
mU = 0.125;
mW = 0.1;
n0 = 0.0004;
np = 0.0375;
nm = 0.5;
gY = 1;
mY = 0.125;

% estoqiometric matrix
   %[W Y U C]
v = [1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    -1 0 -1 1;
    1 0 1 -1;
    -1 0 0 0;
    0 -1 0 0;
    0 0 -1 0;
    0 0 0 -1;
    -1 0 0 0;
    0 -1 0 0;
    0 0 -1 0;
    0 0 0 -1;
    1 0 0 -1;
    0 0 1 -1];


j1 = 1;    % contador para almacenar el numero de reacciones
while t < total_time
    
    % Perturbacion
    if (t > total_time/2)
       mY = 0.5;
    end
    
    % j1p es la reaccion que ocurre (mas o menos) al tiempo de la perturbacion
    if (t > (total_time/2)-10 && t < (total_time/2)+10)
       j1p = j1;
    end

    % Reaction propencities - lamda
    l(1) = mW;            % M synthesis
    l(2) = mY*W;          % Y synthesis 
    l(3) = mU*Y;          % U synthesis
    l(4) = np*U*W;        % complex C formation
    l(5) = n0*C;          % complex C spontaneous unbinding rate
    l(6) = gW*W;          % W degradation 
    l(7) = gY*Y;          % Y degradation
    l(8) = gU*U;          % U degradation
    l(9) = nm*C;          % C degradation
    l(10) = g*W;          % W global dilution
    l(11) = g*Y;          % Y global dilution
    l(12) = g*U;          % U global dilution
    l(13) = g*C;          % C global dilution
    l(14) = gU*C;         % U degradates in the complex form 
    l(15) = gW*C;         % W degradates in the complex form 

    % Calculate total reaction rate
    total_rate = sum(l);
    
    if total_rate == 0
        break; % No more reactions can occur
    end

    % Calculate time until next reaction
    r1 = rand();
    tau = -log(1-r1)/total_rate;

    % Determine which reaction occurs
    r2 = rand();
    l_rates = l/total_rate;  % vector con cada rate para cada reaccion
    reaction = 1;            % elejir una reaccion
    sumR = l_rates(1);       % guardar la suma

    while (sumR < r2)        
        reaction = reaction + 1;
        sumR = sumR + l_rates(reaction);
    end
    
    X = X + v(reaction,:);   % se actualiza el numero de moleculas
    W = X(1);
    Y = X(2);
    U = X(3);
    C = X(4);

    % Update time 
    t = t + tau;
    
    j1 = j1+1;
    % Store species values and time
    X_results(j1,:) = X;
    time_points(j1) = t;
end

% Plot the results for each species
tiledlayout(2,1)
nexttile
plot(time_points(1:j1), X_results((1:j1),:), 'LineWidth', 0.25);
legend({"W", "Y", "U", "C"})
nexttile
plot(time_points(1:j1), X_results((1:j1),2), 'LineWidth', 0.25);
legend({"Y"})


tiledlayout(2,2)
nexttile;
histogram(X_results(1:j1p,1), 'FaceAlpha',0.3, 'Normalization', 'probability');
hold on
histogram(X_results(j1p:j1,1), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"W", "Wp"})
nexttile
histogram(X_results(1:j1p,2), 'FaceAlpha',0.3, 'Normalization', 'probability');
hold on
histogram(X_results(j1p:j1,2), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"Y", "Yp"})
nexttile
histogram(X_results(1:j1p,3), 'FaceAlpha',0.3, 'Normalization', 'probability');
hold on
histogram(X_results(j1p:j1,3), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"U", "Up"})
nexttile
histogram(X_results(1:j1p,4), 'FaceAlpha',0.3, 'Normalization', 'probability');
hold on
histogram(X_results(j1p:j1,4), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"C", "Cp"})




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estado estacionario de Y
n = j1p;
Y_simulation = X_results(1:j1p,2);          % simulacion de Y


Y_molecules = unique(Y_simulation).';   % rango de moleculas de Y
Y_counts = histcounts(Y_simulation);    % conteo de las moleculas de Y
Y_prob = Y_counts/n;                    % probabilidad de cada numero de moleculas en la simulacion


for i = 1:length(Y_molecules)
    Y_cumprob(i) = sum(Y_prob(1:i));
end

Y_distribution = [Y_molecules; Y_cumprob];    % distribucion acumulada de Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Simulation of localy analogous system with fixed imput mU = mU*Y_ss

t = 0;
time_points_analog = zeros(1, 10000000);  

% Initial conditions
W = 0;
Y = 0;
U = 0;
C = 0;
X = [W, Y, U, C]; 

X_results_analog = zeros(10000000, 4); 

% Kinetic parameters same as before
mY = 0.125; % reiniciar mY

% Estoqiometric matrix same as before

j2 = 1;
while t < total_time

    %Perturbation
    if (t>total_time/2)
        mY = 0.5;
    end

    
    if (t > (total_time/2)-10 && t < (total_time/2)+10)
       j2p = j2;
    end

    % Un dato aleatoreo de la simulacion anterior
    % rY = randi(length(Y_simulation));
    % Y_ss = Y_simulation(rY);


    % Usar el vector de la simulacion anterior
    if(j2 < j1p)
         Y_ss = Y_simulation(j2);
    else
         break     
         % si el numero de reacciones supera el numero de reacciones de la simulacion 
         % anterior antes de la perturbacion (j1p) se rompe, aqui no suele
         % pasar 
    end


    % Dato aleatoreo de la distribucion
    % rY = rand();
    % i = 1;
    % while (Y_distribution(2,i) < rY)
    %     i = i +1;
    % end
    % Y_ss = Y_distribution(1,i);
    % % 


    % Reaction propencities - lamda
    l(1) = mW;            % M synthesis
    l(2) = mY*W;          % Y synthesis 

    l(3) = mU*Y_ss;          % U synthesis

    l(4) = np*U*W;        % complex C formation
    l(5) = n0*C;          % complex C spontaneous unbinding rate
    l(6) = gW*W;          % W degradation 
    l(7) = gY*Y;          % Y degradation
    l(8) = gU*U;          % U degradation
    l(9) = nm*C;          % C degradation
    l(10) = g*W;          % W global dilution
    l(11) = g*Y;          % Y global dilution
    l(12) = g*U;          % U global dilution
    l(13) = g*C;          % C global dilution
    l(14) = gU*C;         % like U degradates in the complex form
    l(15) = gW*C;         % like W degradates in the complex form


    % Calculate total reaction rate
    total_rate = sum(l);

    if total_rate == 0
        break; % No more reactions can occur
    end

    % Calculate time until next reaction
    r1 = rand();
    tau = -log(1-r1)/total_rate;

    % Determine which reaction occurs
    r2 = rand();
    l_rates = l/total_rate;
    sumR = l_rates(1);
    reaction = 1;
    while (sumR < r2)
        reaction = reaction + 1;
        sumR = sumR + l_rates(reaction);
    end

    X = X + v(reaction,:);
    W = X(1);
    Y = X(2);
    U = X(3);
    C = X(4);

    % Update time and species counts
    t = t + tau;


    j2 = j2+1;
    % Store species values and time
    X_results_analog(j2,:) = X;
    time_points_analog(j2) = t;
end




% Plot the results for each species
tiledlayout(2,1)
nexttile
plot(time_points_analog(1:j2), X_results_analog((1:j2),:), 'LineWidth', 0.25);
xlim([0 total_time]);
legend({"W", "Y", "U", "C"})
nexttile
plot(time_points_analog(1:j2), X_results_analog((1:j2),2), 'LineWidth', 0.25);
xlim([0 total_time]);
legend({"Y"})


tiledlayout(2,2)
nexttile;   
histogram(X_results_analog(1:j2p,1), 'FaceAlpha',0.3, 'Normalization', 'probability');
hold on
histogram(X_results_analog(j2p:j2,1), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"W", "Wp"})
nexttile
histogram(X_results_analog(1:j2p,2), 'FaceAlpha',0.3, 'Normalization', 'probability');
hold on
histogram(X_results_analog(j2p:j2,2), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"Y", "Yp"})
nexttile
histogram(X_results_analog(1:j2p,3), 'FaceAlpha',0.3, 'Normalization', 'probability', NumBins=max(X_results_analog(1:j2p,3)));
hold on
histogram(X_results_analog(j2p:j2,3), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"U", "Up"})
nexttile
histogram(X_results_analog(1:j2p,4), 'FaceAlpha',0.3, 'Normalization', 'probability');
hold on
histogram(X_results_analog(j2p:j2,4), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"C", "Cp"})
