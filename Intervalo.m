function [ min, max ] = Intervalo( pdf, ToF, porcentaje, tpo_predic, delta_T )

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                           Entradas
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

        %pdf: recibe los datos de la pdf
        %ToF: recibe el tiempo de falla
        %porcentaje: recibe el intervalo de confianza [%]
        %tpo_predic: instante en el que comienza la predicción
        %delta_T: Delta tiempo (entre mediciones)
       
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                           Salidas
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

        %min: límite inferior del intervalo de confianza
        %max: límite superior del intervalo de confianza
        
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                   Iniciación de variables
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

ToF=ToF-tpo_predic;

if ToF < 0
    
    min=0;
    max=0;

else

largo=length(pdf);                                                          %largo del vector de datos.

ToF=round(ToF);                                                             %se necesita el valor en números enteros, por lo que se realiza un redondeo simple

prob=pdf(ToF);                                                              %se inicializa la probabilidad con la probabilidad de estar en el ToF

max=ToF;                                                                    %ambos límites comienzan en la posición del ToF
min=ToF;

%Condiciones limites del intervalo

limite_inf=1;
limite_sup=largo;

for i=1:1:largo
    
    if pdf(i)~= 0
        limite_inf=i;
        break
    end
end

for i=largo:-1:1
    
    if pdf(i)~=0
        limite_sup=i;
        break
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%              Obtención del Intervalo de Confianza
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

for i=1:1:largo
    
    if (prob < porcentaje)
        
        if(min>limite_inf && max<limite_sup)
            
            min=min-1;
            max=max+1;
            
            prob=prob+pdf(min)+pdf(max);
            
        else
            if(min>limite_inf && max>=limite_sup)
                
                min=min-1;
                prob=prob+pdf(min);
            else
                if (min<=limite_inf && max<limite_sup)
                    
                    max=max+1;
                    prob=prob+pdf(max);
                else
                    break;
                end
            end
        end
    else 
        
        break;
    end
end

min=(min+tpo_predic)*delta_T;
max=(max+tpo_predic)*delta_T;

end
end
