function [ x_out ] = regularizacion(xReg,H,kernel,kernel_cdf)

    [N,n_estados]=size(xReg);
 
    %desv = std(xReg);
    %En promedio es mas rapido calcularlo directamente. (Dividido por N-1
    %para que sea insesgado (tal como lo hace matlab).
    desv = sqrt((sum(xReg.^2)-N*mean(xReg).^2)./(N-1));
       
    aleatorio=rand(N,1);        % se crea un valor aleatorio entre 0 y 1
    indice = zeros(N,1);
    
    
    %x_out2 = zeros(N,n_estados);
    %indice_fast = zeros(N,1);
    %[sorted,iSorted] = sort(aleatorio);
    %indice_fast(1) = find(kernel_cdf>=sorted(1),1); 
    %for j=2:N
    %    indice_fast(j) = find(kernel_cdf(indice_fast(j-1):end)>=sorted(j),1)+indice_fast(j-1)-1;   %se buscan los indices de los elementos con valor mayor o igual al n° aleatorio.
    %end
    %x_out(iSorted)=xReg(iSorted)+H.*desv.*kernel(indice_fast,1);
    %indx = sortrows([iSorted indice_fast]);
   
  
    for j=1:N
        indice(j) = find(kernel_cdf>=aleatorio(j),1);   %se buscan los indices de los elementos con valor mayor o igual al n° aleatorio.
        %x_out(j)=xReg(j)+H.*desv.*kernel(indice,1);
    end
    
    x_out=xReg+H.*desv.*kernel(indice,1);
    
end
