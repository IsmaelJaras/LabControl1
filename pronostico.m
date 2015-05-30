 function [prog, s_prog] = pronostico(Ipron,s_prog,mod)
%%
    len_predic = s_prog.h_predic;
    realizaciones = s_prog.realizaciones;
    
    %Estructura que contiene los resultados de cada realizacion
    prog.reali = struct([]);
    prog.mean  = struct('Info','Estructura con resultado promedio de las realizaciones de pronostico');
    prog.mean.soc = zeros(len_predic,1);
    prog.mean.dens = zeros(len_predic,1);
    
    null_SOC_count = zeros(realizaciones,1); %contador de las realizaciones en que no se logró cruzar el umbral de falla

    for i=2:length(Ipron)
        if(Ipron(i)==0)
            Ipron(i) = Ipron(i-1);
        end
    end
    
    lim_sup_SOC = 0.055;
    lim_inf_SOC = 0.045;
    hz_SOC = lim_sup_SOC-lim_inf_SOC;    
    
    %Resampling antes de empezar pronostico con Regularizacion
    [s_prog.part, s_prog.pesos] = resample_sys_dp(s_prog.part, s_prog.pesos);      
    s_prog.obs(:,1) = exp(-(1-s_prog.part(:,2))*mod.g)*(mod.Vo-mod.VL) + mod.VL*(ones(s_prog.npart,1) + mod.alfa*(s_prog.part(:,2) - ones(s_prog.npart,1)) + (1-mod.alfa)*(exp(-mod.beta)-exp(-mod.beta*sqrt( abs(s_prog.part(:,2))))))- Ipron(end)*s_prog.part(:,1);
    s_prog.salida = s_prog.pesos'*s_prog.obs(:,1);
  
    for j=1:s_prog.realizaciones
        %display(sprintf('Realizacion n° %i de etapa de pronostico',j))
 
        prog.reali(j).Isint = markov_I(Ipron,s_prog.h_predic);
        %prog.reali(j).Isint = mean(Ipron)*ones(s_prog.h_predic,1);
        %Isint = I(fp.tpo_predic+1:end);

        % se realiza la predicción a "n" pasos
        [SOC_predic_N, prog.reali(j).x1, prog.reali(j).x2, prog.reali(j).voltage, prog.reali(j).soc] = prediccion(s_prog, mod, prog.reali(j).Isint, len_predic);
        
        %prog.reali(j).soc=(fp.pesos'*SOC_predic_N)';       %SOC esperado en cada instante de tiempo futuro




%--------------------------------------------------------------------------
%                  Generar distribución acumulada
%--------------------------------------------------------------------------
            dist_acum_SOC=zeros(1,len_predic);

%             for jd=1:fp.npart;                                                            %N: Número de Partículas
%                 for i=1:len_predic                                                               %largo: Instantes de Predicción
%                     if (SOC_predic_N(jd,i)<=lim_sup_SOC)
%                         dist_acum_SOC(i)=dist_acum_SOC(i)+fp.pesos(jd)*min(abs(SOC_predic_N(jd,i)-lim_sup_SOC)/hz_SOC,1);
%                     end
%                 end
%             end

            %Cambien el FOR anterior por esto. Anda mas rapido
            for jd=1:s_prog.npart;                                                            %N: Número de Partículas
                [~,b] = find(SOC_predic_N(jd,:)<=lim_sup_SOC);
                dist_acum_SOC(b) = dist_acum_SOC(b)+s_prog.pesos(jd)*min(abs(SOC_predic_N(jd,b)-lim_sup_SOC)/hz_SOC,1);
            end

            
            
            prog.reali(j).CDF = dist_acum_SOC;
%--------------------------------------------------------------------------
%              Generar la densidad de probabilidad
%--------------------------------------------------------------------------
            densidad_SOC=PDF(dist_acum_SOC);

            if sum(densidad_SOC)~=0
                densidad_SOC=densidad_SOC/sum(densidad_SOC);                        % esto es para no tener divisiones por cero.    
            end

            densidad_SOC(1)=0;
            densidad_SOC(len_predic)=0;

            indice_SOC=(1:len_predic)'+s_prog.tpo_predic;
            
            prog.reali(j).tof = indice_SOC'*densidad_SOC';           
            [prog.reali(j).min_IC, prog.reali(j).max_IC ] = Intervalo( densidad_SOC, prog.reali(j).tof, s_prog.confianza, s_prog.tpo_predic, mod.dt );   
      
            prog.reali(j).densidad = densidad_SOC';
            prog.reali(j).nula = 0;
             
            if prog.reali(j).tof==0
                null_SOC_count(j,1)=1;       % revisa si se cruzó o no el umbral de falla
                prog.reali(j).nula = 1;
            end
    end

    contador_null_SOC_total=sum(null_SOC_count(:,1)); %se obtienen cuántas realizaciones no cruzaron el umbral de falla

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%  Promediar los resultados obtenidos para el instante de predicción
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Cuando no se entra a la Hazard Zone.. no se consideran los resultados de esa realización.

%%
%--------------------------------------------------------------------------
%                              SOC
%--------------------------------------------------------------------------

    if realizaciones==contador_null_SOC_total                                    
        total_SOC=1;
    else
        total_SOC=realizaciones-contador_null_SOC_total;
    end

    if realizaciones~=1   % al no cruzar el umbral de falla tanto el valor del ToF como los limites del intervalo de confianza, por lo que se promedian los resultados de las realizaciones que si lograron entrar a la Hazard zone

        prog.mean.tof = sum([prog.reali.tof])/total_SOC;
        prog.mean.min_IC = sum([prog.reali.min_IC])/total_SOC;
        prog.mean.max_IC = sum([prog.reali.max_IC])/total_SOC;

        if contador_null_SOC_total ~= 0  % se suman las predicciones de las realizaciones "válidas"
            for j=1:realizaciones
                if(null_SOC_count(j,1)~=1)
                    prog.mean.soc=prog.mean.soc+prog.reali(j).soc;
                    prog.mean.dens=prog.mean.dens+prog.reali(j).densidad;
                end
            end
            % se promedian las predicciones
            prog.mean.soc=prog.mean.soc/total_SOC;
            prog.mean.dens=prog.mean.dens/total_SOC;
        else
            %Como son estructuras no es llegar y hacer sum().
            %Al hacer [struct.algo] pone en la misma fila todos los
            %elementos de struct(i).algo de i=1:final. Entonces se hace el
            %reshape
            prog.mean.soc =  sum([prog.reali.soc],2)/realizaciones;
            prog.mean.dens = sum([prog.reali.densidad],2)/realizaciones;
        end

    else
        prog.mean.soc = prog.reali.soc;
        prog.mean.tof = prog.reali.tof;
        prog.mean.min_IC = prog.reali.min_IC;
        prog.mean.max_IC = prog.reali.max_IC;
        prog.mean.dens = prog.reali.densidad;
    end  
    
    prog.mean.contNull = null_SOC_count;
end
