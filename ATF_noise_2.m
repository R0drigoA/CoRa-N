clear; clc;
%rng(1)

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


j1 = 1;    % contador para almacenar el numero de ciclos

while t < total_time
    
    if (t > total_time/2)
       mY = 0.5;
    end

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
    l_rates = l/total_rate;  
    reaction = 1;            
    sumR = l_rates(1);       

    while (sumR < r2)        
        reaction = reaction + 1;
        sumR = sumR + l_rates(reaction);
    end
    
    X = X + v(reaction,:);   % la ultima reaccion del ciclo se suma a X y se actualiza el numero de moleculas
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estado estacionario de W

n=j1p;
W_simulation = X_results(1:j1p,1);          % simulacion de Y


W_molecules = unique(W_simulation).';   % rango de moleculas de Y
W_counts = histcounts(W_simulation);    % conteo de las moleculas de Y
W_prob = W_counts/n;                    % probabilidad de cada numero de moleculas en la simulacion


for i = 1:length(W_molecules)
    W_cumprob(i) = sum(W_prob(1:i));
end

W_distribution = [W_molecules; W_cumprob];    % distribucion acumulada de Y


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%% Simulation of localy analogous system without feedback
total_time = 50000;

t = 0;
time_points_analog_v = zeros(1, 10000000);

% Initial conditions
W = 0;
Y = 0;
V = 0;
U = 0;
C = 0;
X = [W, Y, V, U, C]; 
X_results_analog_v = zeros(10000000, 5);

% Kinetic parameters
mY = 0.125;
gV = 1;
mV = 0.125;

% % estoqiometric matrix
v_analog = [1 0 0 0 0;
    0 1 0 0 0;
    0 0 0 1 0;
    -1 0 0 -1 1;
    1 0 0 1 -1;
    -1 0 0 0 0;
    0 -1 0 0 0;
    0 0 -1 0 0;
    0 0 0 -1 0;
    0 0 0 0 -1
    -1 0 0 0 0;
    0 -1 0 0 0;
    0 0 -1 0 0;
    0 0 0 -1 0;
    0 0 0 0 -1;
    1 0 0 0 -1;
    0 0 0 1 -1;
    0 0 1 0 0];


j3 = 1;
j = 1;
while t < total_time
    %Perturbation
    if (t>total_time/2)
        mY = 0.5;
    end

    if (t > (total_time/2)-10 && t < (total_time/2)+10)
       j3p = j3;
    end

    % Obtener una W aleatorea de la distribucion
    % rW = rand();
    % i = 1;
    % while (W_distribution(2,i) < rW)
    %     i = i +1;
    % end
    % W_ss = W_distribution(1,i);  

    % seleccionar una W aleatorea de la simulacion anterior
    rW = randi(j1p);
    W_ss = W_simulation(rW);


    % 
    % if(j3 < j1p)
    %     W_ss = W_simulation(j);
    % else
    %     break 
    % % Se necesita primero hacer que la simulacon anterior sea mucho mas
    % % grande que esta para no romper el ciclo
    % end


    % Reaction propencities - lamda
    l(1) = mW;            % M synthesis
    l(2) = mY*W;          % Y synthesis 
    l(3) = mU*V;          % U synthesis
    l(4) = np*U*W;        % complex C formation
    l(5) = n0*C;          % complex C spontaneous unbinding rate
    l(6) = gW*W;          % W degradation 
    l(7) = gY*Y;          % Y degradation
    l(8) = gV*V;          % V degradation
    l(9) = gU*U;          % U degradation
    l(10) = nm*C;         % C degradation
    l(11) = g*W;          % W global dilution
    l(12) = g*Y;          % Y global dilution
    l(13) = g*V;          % V global dilution
    l(14) = g*U;          % U global dilution
    l(15) = g*C;          % C global dilution
    l(16) = gU*C;         % like U degradates in the complex form
    l(17) = gW*C;         % like W degradates in the complex form
    l(18) = mV*W_ss;


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

    X = X + v_analog(reaction,:);
    W = X(1);
    Y = X(2); 
    V = X(3);
    U = X(4);
    C = X(5);

    % Update time and species counts
    t = t + tau;


    j3 = j3+1;
    % Store species values and time
    X_results_analog_v(j3,:) = X;
    time_points_analog_v(j3) = t;
end



% 
% % Plot the results for each species
tiledlayout(2,1)
nexttile
plot(time_points_analog_v(1:j3), X_results_analog_v((1:j3),:), 'LineWidth', 0.25);
xlim([0 total_time]);
legend({"W", "Y", "V", "U", "C"})
nexttile
plot(time_points_analog_v(1:j3), X_results_analog_v((1:j3),2), 'LineWidth', 0.25);
xlim([0 total_time]);
legend({"Y"})


tiledlayout(2,3)
nexttile;
histogram(X_results_analog_v(1:j3p,1), 'FaceAlpha',0.3, 'Normalization', 'probability');
hold on
histogram(X_results_analog_v(j3p:j3,1), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"W", "Wp"})
nexttile
histogram(X_results_analog_v(1:j3p,2), 'FaceAlpha',0.3, 'Normalization', 'probability');
hold on
histogram(X_results_analog_v(j3p:j3,2), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"Y", "Yp"})
nexttile
histogram(X_results_analog_v(1:j3p,3), 'FaceAlpha',0.3, 'Normalization', 'probability');
hold on
histogram(X_results_analog_v(j3p:j3,3), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"V", "Vp"})
nexttile
histogram(X_results_analog_v(1:j3p,4), 'FaceAlpha',0.3, 'Normalization', 'probability');
hold on
histogram(X_results_analog_v(j3p:j3,4), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"U", "Up"})
nexttile
histogram(X_results_analog_v(1:j3p,5), 'FaceAlpha',0.3, 'Normalization', 'probability');
hold on
histogram(X_results_analog_v(j3p:j3,5), 'FaceAlpha',0.3, 'Normalization', 'probability');
legend({"C", "Cp"})

