close all; clear all;

s_est = struct;
mod = struct;


datos = 2;

data_offset = 0;%Desde que dato se toma con I(1+offset:end)
switch datos 
    case 1
        run datos_bat1_robot
    case 2
        run datos_E1_FUDS
end


%====================================================
%         Inicialización del filtro de partículas
%====================================================

s_est.nestados = 2; %numero de estados
s_est.npart = 40; %numero de particulas
s_est.CI = zeros(1,s_est.nestados); %numero de particulas
s_est.part = zeros(s_est.npart,s_est.nestados); %particulas 
s_est.obs = zeros(s_est.npart,1); %particulas
s_est.pesos = ones(s_est.npart,1)./s_est.npart; %pesos
s_est.salida = 4;

%s_est.tpo_predic = floor(l_datos*0.65);
s_est.tpo_predic = 4100;
s_est.h_predic = length(I)-s_est.tpo_predic ;
s_est.confianza = 0.95;
s_est.realizaciones = 25;

%====================================================
%               Inicializacion del f
%====================================================

mod.E0 = 0.85;
s_est.CI(1) = r0;

s_est.part(:,1) =  s_est.CI(1)+0.005.*randn(s_est.npart,1);
s_est.part(:,2) = mod.E0 + 0.085*rand(s_est.npart,1);

s_est.pesos = ones(s_est.npart,1)./s_est.npart;

%pesos=ones(l_datos,fp.npart)/fp.npart;


%====================================================
%               SOC
%====================================================

soc_counting     = 1-cumsum(V.*I)./sum(V.*I);
soc_real = soc_real(1+data_offset:end);
%%
%====================================================
%               Estimación
%====================================================

%Offset: desde donde parten los datos, el offset del principio estaba 
%asignado por Daniel pero no me funcionó
offset = 0;
V = V(1+offset:end);
I = I(1+offset:end);
soc_real = soc_real(1+offset:end);
soc_counting = soc_counting(1+offset:end);

%Función de estimación (cambiar gráficos según método):
%FP base: estimacion2v2(s_est,mod,V,I)
%FP+OFCL: estimacion2v2_ofcl(s_est,mod,V,I)
[s_est, soc_filtrado, v_estim, imp(1:s_est.tpo_predic), s_est_hist] = estimacion2v2_ofcl(s_est,mod,V,I);


%resampling de partículas históricas para graficar
for i = 1:s_est.tpo_predic
    [xi, ~] = resample_sys_dp(squeeze(s_est_hist.part(i,:,:)), s_est_hist.pesos(i,:)');
    s_est_hist.part(i,:,:) = xi;
end

%Gráfico de estimación cuando no se desea ejecutar el pronóstico.
%run graficos_estimacion_pf_ofcl

disp('Media:')
disp(mean(squeeze(s_est_hist.part(s_est.tpo_predic,:,:))))
disp('Covarianza:')
disp(cov(squeeze(s_est_hist.part(s_est.tpo_predic,:,:))))
disp('Error SOC:')
disp(soc_counting(s_est.tpo_predic)-mean(squeeze(s_est_hist.part(s_est.tpo_predic,:,2))))


%%{ 
%%
%====================================================
%               Pronóstico
%====================================================

Ipron = I(1:s_est.tpo_predic);

s_prog = s_est;

%fprog: estructura usada para filtro de particulas en pronostico
%prog: los resultados
[s_result, s_prog] = pronostico(Ipron,s_prog,mod);   



%%
%Resultados
jit5 = find(cumsum(s_result.mean.dens)>0.05,1) + s_est.tpo_predic;
jit15 = find(cumsum(s_result.mean.dens)>0.15,1) + s_est.tpo_predic;

display('===========================')
display(' ')
display(sprintf('EOD >> %4.0f',round(s_result.mean.tof)))
display(sprintf('IC >> [%4.0f ; %4.0f]',round(s_result.mean.min_IC), round(s_result.mean.max_IC)))
display(sprintf('JIT_5 >> %4.0f',round(jit5)))
display(sprintf('JIT_15 >> %4.0f',round(jit15)))

display(' ')
display('===========================')



%% Graficos Estimacion y Pronostico
graficar = 4;

%REALIZACION A GRAFICAR!
reali_plot = 3;
%================%

close all
switch graficar
    case 1
        run graficos_color.m
    case 2
        run graficos_bn.m %Blanco y negro para el paper
    case 3
        run graficos_color_pf.m
    case 4
        run graficos_color_pf_ofcl.m
    otherwise
end

annotation(figure(3),'textbox',...
    [0.039 0.861 0.039 0.099],...
    'String',{'c)'},...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'LineStyle','none');

annotation(figure(2),'textbox',...
    [0.039 0.861 0.039 0.099],...
    'String',{'a)'},...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'LineStyle','none');

annotation(figure(2),'textbox',...
    [0.039 0.389 0.039 0.099],...
    'String',{'b)'},...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'LineStyle','none');

annotation(figure(4),'textbox',...
    [0.039 0.861 0.039 0.099],...
    'String',{'d)'},...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'LineStyle','none');

%}
