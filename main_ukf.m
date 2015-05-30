%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main modificado de la versi�n de Daniel Pola para hacer uso del Filtro de
% Kalman Unscented.
% Por Carlos Tampier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%      Inicializaci�n de la Estructura de Datos
%====================================================
% Variables para estimaci�n. Condiciones iniciales y part�culas para
% predicci�n

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

% Variables UKF

s_est.media = zeros(s_est.nestados,1);
s_est.cov = zeros(s_est.nestados,s_est.nestados);


%====================================================
%               Condiciones Iniciales
%====================================================

mod.E0 = 0.85;%condici�n inicial de SOC
s_est.CI(1) = r0;%condici�n inicial de Zin

%====================================================
%               SOC
%====================================================

soc_counting     = 1-cumsum(V.*I)./sum(V.*I);
soc_real = soc_real(1+data_offset:end);
%%
%====================================================
%               Estimaci�n
%====================================================

%Offset: desde donde parten los datos, el offset del principio estaba 
%asignado por Daniel pero no me funcion�
offset = 0; 
V = V(1+offset:end);
I = I(1+offset:end);
soc_real = soc_real(1+offset:end);
soc_counting = soc_counting(1+offset:end);

%Funci�n de estimaci�n (cambiar gr�ficos seg�n m�todo):
%UKF: estimacion_ukf(s_est,mod,V,I)
%UKF+OFCL: estimacion_ukf_ofcl5(s_est,mod,V,I);
%[s_est, soc_filtrado, v_estim, imp(1:s_est.tpo_predic), s_est_hist] = estimacion_ukf_ofcl5(s_est,mod,V,I);
[s_est, soc_filtrado, v_estim, imp(1:s_est.tpo_predic), s_est_hist] = estimacion_ukf(s_est,mod,V,I);

%desviaci�n est�ndar de los estados en el tiempo
std_zin = nan(s_est.tpo_predic,1);
std_soc = nan(s_est.tpo_predic,1);
for i = 1:s_est.tpo_predic
    cov = squeeze(s_est_hist.cov(i,:,:));
    std_zin(i) = sqrt(cov(1,1));
    std_soc(i) = sqrt(cov(2,2));
end

%Gr�fico de estimaci�n cuando no se desea ejecutar el pron�stico.
%run graficos_estimacion_ukf_ofcl

disp('Media:')
disp(s_est.media')
disp('Covarianza:')
disp(s_est.cov)
disp('Error SOC:')
disp(soc_counting(s_est.tpo_predic)-s_est.media(2))

media_ukf = s_est.media;
cov_ukf = s_est.cov;

%save('perf_ukf_085_4100.mat','media_ukf','cov_ukf')

%%{
%%
%====================================================
%        Generaci�n de Part�culas
%====================================================

s_est.part = max(mvnrnd(repmat(s_est.media',s_est.npart,1),s_est.cov),0);
s_est.pesos = ones(s_est.npart,1)./s_est.npart;


%%
%====================================================
%               Pron�stico
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
reali_plot = 2;
%================%

close all
switch graficar
    case 1
        run graficos_color.m
    case 2
        run graficos_bn.m %Blanco y negro para el paper
    case 3
        run graficos_color_ukf.m %Graficos a color para UKF
    case 4
        run graficos_color_ukf_ofcl.m %Graficos a color para UKF+OFCL
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
