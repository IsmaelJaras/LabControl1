function [ pdf ] = PDF( cdf )

pdf=diff(cdf);

largo=length(pdf);

aux=0;  

%--------------------------------------------------------------------------
%    CONDICIÓN DE ESPECIAL: CDF=0 EN TODO INSTANTE (NO CRUZA UMBRAL)
%--------------------------------------------------------------------------

if sum(pdf)~=0
    
%--------------------------------------------------------------------------
%       CONDICIÓN DE BORDE: ÚLTIMO TÉRMINO NEGATIVO
%--------------------------------------------------------------------------

    if pdf(largo) < 0 

        pdf(largo)=0;

    end

%--------------------------------------------------------------------------
%       CONDICIÓN DE BORDE: PRIMER TÉRMINO NEGATIVO
%--------------------------------------------------------------------------

    if pdf(1) < 0

        for j=2:1:largo

                if pdf(j)>=0              %Encuentro la posición del siguiente valor positivo

                    aux=j;
                    break
                end
        end

           % y= m*x + n
           % m = (y2-y1)/(x2-x1)
           % n = y1-m*x1 || n= y2-m*x2
           %en este caso:
           % y1=pdf(i) ; y2 = pdf(aux)
           % x1=i ; x2 = aux

           m = (pdf(aux)-0)/(aux-0);
           n = pdf(aux)-m*aux;

           % Hay que reemplazar los valores negativos.

            for j=1:1:aux-1

           pdf(j) = m*j + n;

            end 

    end  

%--------------------------------------------------------------------------
%               TÉRMINOS FUERA DE LAS CONDICIONES DE BORDE
%--------------------------------------------------------------------------

    for i=1:1:largo-1

        if pdf(i+1) < 0                   %Encuentro la posición del valor positivo(o cero) anterior al valor negativo

            for j=i+1:1:largo

                if pdf(j)>=0              %Encuentro la posición del siguiente valor positivo

                    aux=j;
                    break
                end
            end


            %hay que encontrar la ecuación de la recta entre los dos valores positivos

                        m = (pdf(aux)-pdf(i))/(aux-i);
                        n = pdf(i)-m*i;

            % Hay que reemplazar los valores negativos.

            for j=i+1:1:aux-1

                pdf(j) = m*j + n;

            end

        end
    end
end
end
