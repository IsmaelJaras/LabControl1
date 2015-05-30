function [ SOC_predic_N, x1, x2, v_estp, soc_estp ] = prediccion(fp, mod, Isint,len_predic)
 

%(x1,x2,I_pre,delta_T,largo,valores_filtro,E_crit)

SOC_predic_N=zeros(fp.npart,len_predic);    %Vector que contiene el valor del SOC por cada particula en cada instante
x1 = zeros(fp.npart,len_predic);
x2 = zeros(fp.npart,len_predic);
obs = zeros(fp.npart,len_predic);
v_estp = zeros(1,len_predic);

%--------------------------------------------------------------------------
%                         Kernel y H Regularizacion
%--------------------------------------------------------------------------
%Estos calculos estan fuera de la funcion de regularizacion por que al
%tener len_predic*n_realizaciones => 4500*40 = 180.000 el mismo calculo que
%toma un buen tiempo
    C_n_estados = 1;
    A=(8.*(C_n_estados)^(-1).*(fp.nestados+4).*(2.*sqrt(pi)).^fp.nestados).^(1./(fp.nestados+4));
    H=A.*fp.npart.^(-1./(fp.nestados+4));
    %%%%ajuste%%%%%%%%%%%%%
    H=H/30; %13
    %%%%%%%%%%%%%%%%%%%%%%%
    %era kernel=[((1:201)'-1)*0.01-1,sqrt(1-((-1:0.01:1)'.^2))];
    kernel=[(-1:0.01:1)' , sqrt(1-((-1:0.01:1)'.^2))];
    kernel_cdf = cumsum(kernel(:,2))/max(sum(kernel(:,2)));


%--------------------------------------------------------------------------
%                         Predicción
%--------------------------------------------------------------------------
    x1(:,1) = fp.part(:,1);
    x2(:,1) = fp.part(:,2);

    t=1;
    x1(:,t)=regularizacion(x1(:,t),H,kernel,kernel_cdf);
    x2(:,t)=regularizacion(x2(:,t),H,kernel,kernel_cdf);
    

    obs(:,t) = mod.VL*(1 + mod.alfa*(x2(:,t) - 1) + (1-mod.alfa)*(exp(-mod.beta)-exp(-mod.beta*sqrt( abs(x2(:,t))))))- Isint(t)*x1(:,t);
    x2(:,t+1) = abs(x2(:,t) - obs(:,t)*Isint(t)*mod.dt/mod.Ecrit);
    x1(:,t+1) = x1(:,t);
    
    v_estp(t) = sum(obs(:,t))/fp.npart;
    SOC_predic_N(:,2)=x2(:,t);
    
    %separe el t=1 del ciclo
    for t=2:len_predic-1
       
        x1(:,t) = regularizacion(x1(:,t),H,kernel,kernel_cdf);
        x2(:,t) = regularizacion(x2(:,t),H,kernel,kernel_cdf);

        obs(:,t) = mod.VL*(1 + mod.alfa*(x2(:,t) - 1) + (1-mod.alfa)*(exp(-mod.beta)-exp(-mod.beta*sqrt( abs(x2(:,t))))))- Isint(t)*fp.part(:,1);
        
        x1(:,t+1) = x1(:,t);
        x2(:,t+1) = abs(x2(:,t) - obs(:,t)*Isint(t)*mod.dt/mod.Ecrit);

        v_estp(t) = sum(obs(:,t))/fp.npart;
        %Deberia ser fp.pesos'*obs(:,t) pero con regularizacion hay un
        %resampling al principio y todos los pesos son iguales

    % si cruza el valor cero, la partícula se queda con un valor constante para el SOC  
        for n=1:fp.npart
            if SOC_predic_N(n,t)<=0.005 %&& t>2 Al quitar este t>2 puse el caso inicial al principio
                x2(n,t)=0.01;
            end
            SOC_predic_N(n,t+1)=x2(n,t);
        end
    end
    
    soc_estp = (fp.pesos'*SOC_predic_N)';       %SOC esperado en cada instante de tiempo futuro
end
