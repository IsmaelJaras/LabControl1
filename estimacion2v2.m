function [ s_est, soc, vout_estim, imp, s_est_hist] = estimacion2v2( s_est, mod, V, I)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementaci�n de FP para estimaci�n del SOC de bater�as. Trabajo en
% base al Filtro de Part�culas de Daniel Pola (S�lo modelo alterado).
% Por Carlos Tampier.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Vfp = V(j); Ifp = I(j);

% g1 = mod.g;
% g2 = mod.g2;
% alfa = mod.alfa;
% beta = mod.beta;
% Vo = mod.VL;
%observacion = exp((x2-1).*g1)./g2 + Vo*(1+alfa*(x2-1)+(1-alfa)*(exp(-beta)-exp(-beta*sqrt(x2)))) - I*x1;

% Inicializaci�n de vectores y estructuras
vout_estim = 4*ones(s_est.tpo_predic,1);
imp = s_est.CI(1)*ones(s_est.tpo_predic,1);
soc = mod.E0*ones(s_est.tpo_predic,1);
s_est_hist = struct;
s_est_hist.part = nan(s_est.tpo_predic,s_est.npart,s_est.nestados);
s_est_hist.pesos = nan(s_est.tpo_predic,s_est.npart);


for h=1:s_est.tpo_predic
    
    %State-Space Transition
    s_est.part(:,2) = s_est.part(:,2) - (exp(-(1-s_est.part(:,2))*mod.g)*(mod.Vo-mod.VL) + mod.VL*(1 + mod.alfa*(s_est.part(:,2) - 1) + ...
        (1-mod.alfa)*(exp(-mod.beta)-exp(-mod.beta*sqrt(s_est.part(:,2) ))))- I(h)*s_est.part(:,1) ).*I(h)*mod.dt./mod.Ecrit;
    s_est.part(:,1) = s_est.part(:,1) + randn(s_est.npart, 1).*mod.std_x(1);
    s_est.part(:,2) = s_est.part(:,2) + randn(s_est.npart, 1).*mod.std_x(2);%Separado es mas rapido

    s_est.part(:,2) = abs(s_est.part(:,2));

    %Measurement Likelihood
    e2 = V(h) - (exp(-(1-s_est.part(:,2))*mod.g)*(mod.Vo-mod.VL) + mod.VL*(1 + mod.alfa*(s_est.part(:,2) - 1) + ...
        (1-mod.alfa)*(exp(-mod.beta)-exp(-mod.beta*sqrt(s_est.part(:,2) ))))- I(h)*s_est.part(:,1) );

    %Actualizacion de los pesos
    s_est.pesos = s_est.pesos.*(1/sqrt(2*pi*mod.std_obs^2)*exp(-(e2.^2)/(2*mod.std_obs^2)));
    s_est.pesos = s_est.pesos./sum(s_est.pesos); % normalizacion de los pesos

    %Resampling
    Neff = 1/sum(s_est.pesos.^2);

    if (Neff<0.80*s_est.npart)
        [s_est.part, s_est.pesos] = resample_sys_dp(s_est.part, s_est.pesos);  
    end

    %No tener particulas negativas
    s_est.part = abs(s_est.part);

    s_est.obs(:,1) = exp((s_est.part(:,2)-1).*mod.g).*(mod.Vo-mod.VL) + mod.VL*(1 + mod.alfa*(s_est.part(:,2)-1) + (1-mod.alfa)*(exp(-mod.beta)-exp(-mod.beta*sqrt( s_est.part(:,2))))) - I(h)*s_est.part(:,1);

    %Salidas de la funcion de estimacion
    vout_estim(h) = s_est.pesos'*s_est.obs(:,1);
    imp(h) = (s_est.pesos'*s_est.part(:,1));
    soc(h) = (s_est.pesos'*s_est.part(:,2));
    
    s_est_hist.part(h,:,:) = s_est.part;
    s_est_hist.pesos(h,:) = s_est.pesos';
    
    %Reduccion de la varianza del ruido
    if(mod.std_x(2) > 2e-4 && h>200)% && soc(h)<0.97
        mod.std_x(2)=mod.std_x(2)/(1.01);
        %display(h);
    end


end

end
