function [ Isint ] = markov_I( Ipron, h_predic)
%%

    n=2;%h_predic=fp.h_predic;
    Is = smooth(Ipron(1:end),25);
    
    ventana_tiempo = 150;
    n_intervalos = length(Is)/ventana_tiempo;
    resto = rem(length(Is), ventana_tiempo);
    
    ti = 1;
    tf = (ti-1) + ventana_tiempo + resto;
    
    Istate = 0;
    Pewma = zeros(n,n);
    alfa = 0.65;
    
    for intervalo = 1:n_intervalos
        
        s = Is(ti:tf) - mean(Is(ti:tf));
    %Cuantizar s

        gran = max(s);
        peq = min(s);

        x = zeros(length(s),1);
        xlvl = zeros(length(s),1);

        niveles = linspace(peq,gran,n);
        nivelesQ = (1:n);

        for i=1:length(s)
            [~,j] = min(abs(s(i)-niveles));
            x(i) = niveles(j);
            xlvl(i) = nivelesQ(j);
        end

        %Q matriz con niveles de corriente y su estado correspondiente
       

    %Probabilidades de Transicion
    
    
    %POR SI ACASO: ESTA MATRIZ P ESTA TRASPUESTA; QS TAMBIEN ESTA
    %TRASPUESTA. SUMAN 1 EN LAS COLUMNAS.
        P = zeros(n,n);

        for i=1:length(s)-1
            P(xlvl(i+1),xlvl(i))=P(xlvl(i+1),xlvl(i))+1;
        end

        for i=1:length(P)
            if(sum(P(:,i))>0)
                %P(i,:) = P(i,:)/sum(P(i,:));
                P(:,i) = P(:,i)/sum(P(:,i));

            end
        end


    %Generar CM
%     Q = [niveles' nivelesQ'];
%         CI = xlvl(1);  

%         Qs = cumsum(P);
% 
%         x = zeros(largo,1);
% 
%         state = CI*ones(largo,1);

%         for i=2:largo
%             x(i) = rand;
%             for j=1:largo
%                 if(x(i)<Qs(j,state(i-1)))
%                     state(i) = j ;       
%                     break
%                 end
%             end
%         end

%         Isint = zeros(largo,1);
%         for jj = 1:largo;
%             Isint(jj) = Q(state(jj),1)+mean(Ipron);
%         end
        
        ti = tf+1;
        tf = tf + ventana_tiempo;
        
        if(intervalo~=1)
            Istate = (1-alfa).*Istate+alfa.*niveles';
            Pewma = (1-alfa).*P+Pewma*alfa;
        else
            Istate = niveles';
            Pewma = P;
        end
        
        Q = [Istate nivelesQ'];
    end
    
    %%
    CI = xlvl(1);
    
    Qs = cumsum(Pewma);
    x = zeros(h_predic,1);
    state = CI*ones(h_predic,1); 
    
    for i=2:h_predic
        x(i) = rand;
        for j=1:length(Qs)
            if(x(i)<Qs(j,state(i-1)))
                state(i) = j ;       
                break
            end
        end
    end
    
    Isint = zeros(h_predic,1);
    for jj = 1:h_predic;
        Isint(jj) = Q(state(jj),1)+mean(Ipron);
    end


end
