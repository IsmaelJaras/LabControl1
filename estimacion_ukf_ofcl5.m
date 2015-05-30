function [ s_est, soc, vout_estim, imp, s_est_hist] = estimacion_ukf_ofcl5( s_est, mod, V, I)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementación de UKF para estimación del SOC de baterías. Trabajo en
% base al Filtro de Partículas de Daniel Pola.
% Por Carlos Tampier.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Inicialización de los datos de salida
vout_estim = 4*ones(s_est.tpo_predic,1);
imp = s_est.CI(1)*ones(s_est.tpo_predic,1);
soc = mod.E0*ones(s_est.tpo_predic,1);
s_est_hist = struct;
s_est_hist.cov = nan(s_est.tpo_predic,2,2);
s_est_hist.ofcl = zeros(s_est.tpo_predic,1);


% Valores iniciales media, matriz de covarianza y observación predicha
x_med = [s_est.CI(1); mod.E0];
Rxx = [0.005^2, 0; 0, 0.05^2];

% Matriz de varianza para ruido de observación 
Rvv = mod.std_obs^2;

% Parámetros del UKF
ukf = struct;
ukf.alpha = 1;
ukf.beta = 0;
ukf.kappa = 1;

ukf.lambda = ukf.alpha^2*(s_est.nestados+ukf.kappa) - s_est.nestados;
ukf.gamma = sqrt(s_est.nestados+ukf.lambda);

ukf.W0m = ukf.lambda/(s_est.nestados + ukf.lambda);
ukf.W0c = ukf.lambda/(s_est.nestados + ukf.lambda) + (1-ukf.alpha^2+ukf.beta);
ukf.W1 = 1/(2*(s_est.nestados + ukf.lambda));

% Umbral de error para OFCL
eth = 0.15;
err_acum = 0;
% Parámetros OFCL
%p_ofcl = [0.99,0.98];
%q_ofcl = [1.1,1.01];
p_ofcl = [0.99,0.98];
q_ofcl = [1.1,1.01];
min_std = [1e-5, 1e-5];
h_min = 5;

% Registro de activación OFCL

% Iteraciones del UKF
for h=1:s_est.tpo_predic
    
    % Cálculo de Sigma Points
    Sxx = chol(Rxx,'lower');
    x_sp = [x_med, x_med + ukf.gamma*Sxx(:,1), x_med + ukf.gamma*Sxx(:,2), x_med - ukf.gamma*Sxx(:,1), x_med - ukf.gamma*Sxx(:,2)];
    % Predicción del Estado
    x_sp(2,:) = x_sp(2,:) - (exp(-(1-x_sp(2,:))*mod.g)*(mod.Vo-mod.VL) + mod.VL*(1 + mod.alfa*(x_sp(2,:) - 1) + ...
        (1-mod.alfa)*(exp(-mod.beta)-exp(-mod.beta*sqrt(x_sp(2,:) ))))- I(h)*x_sp(1,:) ).*I(h)*mod.dt./mod.Ecrit;
    x_med = ukf.W0m*x_sp(:,1) + ukf.W1*x_sp(:,2) + ukf.W1*x_sp(:,3) + ukf.W1*x_sp(:,4) + ukf.W1*x_sp(:,5);
    
    % Matriz de Covarianza para Ruido de Proceso
    Rww = [mod.std_x(1)^2, 0; 0, mod.std_x(2)^2];
    
    % Predicción de la Covarianza
    Rxx = ukf.W0c*(x_sp(:,1)-x_med)*(x_sp(:,1)-x_med)'+Rww;
    for i = 2:2*s_est.nestados+1
        Rxx = Rxx + ukf.W1*(x_sp(:,i)-x_med)*(x_sp(:,i)-x_med)';
    end
    
    % Cálculo de Sigma Points de Observación
    Sxx = chol(Rxx,'lower');
    x_sp = [x_med, [x_med,x_med] + ukf.gamma*Sxx, [x_med,x_med] - ukf.gamma*Sxx];
    y_sp = (exp(-(1-x_sp(2,:))*mod.g)*(mod.Vo-mod.VL) + mod.VL*(1 + mod.alfa*(x_sp(2,:) - 1) + ...
        (1-mod.alfa)*(exp(-mod.beta)-exp(-mod.beta*sqrt(x_sp(2,:) ))))- I(h)*x_sp(1,:) );
    
    % Predicción de Observación
    y_est = ukf.W0m*y_sp(:,1) + ukf.W1*y_sp(:,2) + ukf.W1*y_sp(:,3) + ukf.W1*y_sp(:,4) + ukf.W1*y_sp(:,5);
    
    % Autocovarianza Residual
    Ryy = ukf.W0c*(y_sp(1)-y_est)^2+Rvv;
    for i = 2:2*s_est.nestados+1
        Ryy = Ryy + ukf.W1*(y_sp(i)-y_est)^2;
    end
    %Ryy = sum([ukf.W0c, ukf.W1*ones(1, 2*s_est.nestados)].*(y_sp - y_est).^2) + Rvv;
    
    % Ganancia
    Rxy = ukf.W0c*(x_sp(:,1)-x_med)*(y_sp(1)-y_est);
    for i = 2:2*s_est.nestados+1
        Rxy = Rxy + ukf.W1*(x_sp(:,i)-x_med)*(y_sp(i)-y_est);
    end
    %Rxy = sum([ukf.W0c*ones(2,1), ukf.W1*ones(2, 2*s_est.nestados)].*(x_sp - repmat(x_med,1,2*s_est.nestados + 1)).*repmat(y_sp - y_est,2,1),2);
    K = Rxy/Ryy;
    
    % Actualización
    err_obs = V(h) - y_est;
    x_med = x_med + K*(err_obs);
    Rxx = Rxx - K*Ryy*K';

    %Salidas de la funcion de estimacion
    vout_estim(h) = exp((x_med(2)-1).*mod.g).*(mod.Vo-mod.VL) + mod.VL*(1 + mod.alfa*(x_med(2)-1) + (1-mod.alfa)*(exp(-mod.beta)-exp(-mod.beta*sqrt( x_med(2))))) - I(h)*x_med(1);
    imp(h) = x_med(1);
    soc(h) = x_med(2);
    
    %Registro histórico
    s_est_hist.cov(h,:,:) = Rxx;
    
    
    % OFCL
    if h>h_min
        
        % Ajuste del ruido asociado a la impedancia interna
        err_acum = err_acum + abs(err_obs);
        if  err_acum <= eth
            mod.std_x = max(mod.std_x.*p_ofcl,min_std);
        else
            err_acum = 0;
            s_est_hist.ofcl(h) = 1;
            mod.std_x = mod.std_x.*q_ofcl;
        end
    end

end

s_est.media = x_med;
s_est.cov = Rxx;

end
