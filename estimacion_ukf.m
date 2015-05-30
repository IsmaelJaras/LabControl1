function [ s_est, soc, vout_estim, imp, s_est_hist] = estimacion_ukf( s_est, mod, V, I)

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


% Valores iniciales media, matriz de covarianza y observación predicha
x_med = [s_est.CI(1); mod.E0];
Rxx = [0.005^2, 0; 0, 0.05^2];

% Matriz de varianza para ruido de observación 
Rvv = mod.std_obs^2;
% Matriz de Covarianza para Ruido de Proceso
Rww = [mod.std_x(1)^2, 0; 0, mod.std_x(2)^2];

% Parámetros del UKF
ukf_alpha = 1;
ukf_beta = 0;
ukf_kappa = 1;

ukf_lambda = ukf_alpha^2*(s_est.nestados+ukf_kappa) - s_est.nestados;
ukf_gamma = sqrt(s_est.nestados+ukf_lambda);

W0m = ukf_lambda/(s_est.nestados + ukf_lambda);
W0c = ukf_lambda/(s_est.nestados + ukf_lambda) + (1-ukf_alpha^2+ukf_beta);
W1 = 1/(2*(s_est.nestados + ukf_lambda));

% Iteraciones del UKF
for h=1:s_est.tpo_predic
    %disp(h)
    
    % Cálculo de Sigma Points
    Sxx = chol(Rxx,'lower');
    x_sp = [x_med, x_med + ukf_gamma*Sxx(:,1), x_med + ukf_gamma*Sxx(:,2), x_med - ukf_gamma*Sxx(:,1), x_med - ukf_gamma*Sxx(:,2)];
    % Predicción del Estado
    x_sp(2,:) = x_sp(2,:) - (exp(-(1-x_sp(2,:))*mod.g)*(mod.Vo-mod.VL) + mod.VL*(1 + mod.alfa*(x_sp(2,:) - 1) + ...
        (1-mod.alfa)*(exp(-mod.beta)-exp(-mod.beta*sqrt(x_sp(2,:) ))))- I(h)*x_sp(1,:) ).*I(h)*mod.dt./mod.Ecrit;
    x_med = W0m*x_sp(:,1) + W1*x_sp(:,2) + W1*x_sp(:,3) + W1*x_sp(:,4) + W1*x_sp(:,5);
    
    % Predicción de la Covarianza
    Rxx = W0c*(x_sp(:,1)-x_med)*(x_sp(:,1)-x_med)'+Rww;
    for i = 2:2*s_est.nestados+1
        Rxx = Rxx + W1*(x_sp(:,i)-x_med)*(x_sp(:,i)-x_med)';
    end
    
    % Cálculo de Sigma Points de Observación
    Sxx = chol(Rxx,'lower');
    x_sp = [x_med, [x_med,x_med] + ukf_gamma*Sxx, [x_med,x_med] - ukf_gamma*Sxx];
    y_sp = (exp(-(1-x_sp(2,:))*mod.g)*(mod.Vo-mod.VL) + mod.VL*(1 + mod.alfa*(x_sp(2,:) - 1) + ...
        (1-mod.alfa)*(exp(-mod.beta)-exp(-mod.beta*sqrt(x_sp(2,:) ))))- I(h)*x_sp(1,:) );
    
    % Predicción de Observación
    y_est = W0m*y_sp(:,1) + W1*y_sp(:,2) + W1*y_sp(:,3) + W1*y_sp(:,4) + W1*y_sp(:,5);
    
    % Autocovarianza Residual
    Ryy = W0c*(y_sp(1)-y_est)^2+Rvv;
    for i = 2:2*s_est.nestados+1
        Ryy = Ryy + W1*(y_sp(i)-y_est)^2;
    end
    %Ryy = sum([W0c, W1*ones(1, 2*s_est.nestados)].*(y_sp - y_est).^2) + Rvv;
    
    % Ganancia
    Rxy = W0c*(x_sp(:,1)-x_med)*(y_sp(1)-y_est);
    for i = 2:2*s_est.nestados+1
        Rxy = Rxy + W1*(x_sp(:,i)-x_med)*(y_sp(i)-y_est);
    end
    %Rxy = sum([W0c*ones(2,1), W1*ones(2, 2*s_est.nestados)].*(x_sp - repmat(x_med,1,2*s_est.nestados + 1)).*repmat(y_sp - y_est,2,1),2);
    K = Rxy/Ryy;
    
    % Actualización
    x_med = x_med + K*(V(h) - y_est);
    Rxx = Rxx - K*Ryy*K';

    %Salidas de la funcion de estimacion
    vout_estim(h) = exp((x_med(2)-1).*mod.g).*(mod.Vo-mod.VL) + mod.VL*(1 + mod.alfa*(x_med(2)-1) + (1-mod.alfa)*(exp(-mod.beta)-exp(-mod.beta*sqrt( x_med(2))))) - I(h)*x_med(1);
    imp(h) = x_med(1);
    soc(h) = x_med(2);
    
    %Registro histórico
    s_est_hist.cov(h,:,:) = Rxx;
    

end

s_est.media = x_med;
s_est.cov = Rxx;

end
