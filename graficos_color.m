%=========================%
% Intervalos de Confianza %
%=========================%

lim_inf_x2 = zeros(length(s_result.reali(reali_plot).x2),1);
lim_sup_x2 = lim_inf_x2;
lim_inf95_x2 = lim_inf_x2;
lim_sup95_x2 = lim_inf_x2;

particulas_predic = sort(s_result.reali(reali_plot).x2,1);
impedancia_predic = sort(s_result.reali(reali_plot).x1,1);

impedancia_media = (s_est.pesos'*s_result.reali(reali_plot).x1)';

for j=1:length(lim_inf_x2)
    
    x1Min95 = impedancia_predic(2,j);
    x1Max95 = impedancia_predic(end-1,j);
    x1Min = impedancia_predic(1,j);
    x1Max = impedancia_predic(end,j);

    if(x1Min95>0.03)
        lim_inf_x1(j) = x1Min;
        lim_inf95_x1(j) = x1Min95;
    else
        lim_inf_x1(j) = NaN;
        lim_inf95_x1(j) = NaN;
    end
    if(x1Max95>0.03)
        lim_sup_x1(j) = x1Max;
        lim_sup95_x1(j) = x1Max95;
    else
        lim_sup_x1(j) = NaN;
        lim_sup95_x1(j) = NaN;
    end
    
    x2Min95 = particulas_predic(2,j);
    x2Max95 = particulas_predic(end-1,j);
    x2Min = particulas_predic(1,j);
    x2Max = particulas_predic(end,j);

    if(x2Min95>0.03)
        lim_inf_x2(j) = x2Min;
        lim_inf95_x2(j) = x2Min95;
    else
        lim_inf_x2(j) = NaN;
        lim_inf95_x2(j) = NaN;
    end
    if(x2Max95>0.03)
        lim_sup_x2(j) = x2Max;
        lim_sup95_x2(j) = x2Max95;
    else
        lim_sup_x2(j) = NaN;
        lim_sup95_x2(j) = NaN;
    end
end







%Colores
rojo = [0.89, 0, 0.13];
morado_claro = [0.6, 0.4, 0.8];
naranjo  =  [0.91, 0.41, 0.17];
verde_claro  =  [0.13, 0.55, 0.13];
verde_osc =  [0, 0.27, 0.13];
gris = [0.3804    0.3804    0.3804];
color_densidad = [0.34 0.63 0.83];



%Posicion de las figuras
fig2 = figure(2);
fig3 = figure(3);
fig4 = figure(4);
scnsize = get(0,'ScreenSize');
borders = get(fig2,'OuterPosition') - get(fig2,'Position');
edge = -borders(1)/2;
pos2 = [edge,scnsize(4) * (1/4),scnsize(3)/2 - edge,scnsize(4)*3/4];
pos3 = [scnsize(3)/2 + edge,pos2(2),pos2(3),pos2(4)/2];
pos4 = [scnsize(3)/2 + edge,scnsize(4)/1.6,pos2(3),pos2(4)/2];
set(fig2,'OuterPosition',pos2)
set(fig3,'OuterPosition',pos3)
set(fig4,'OuterPosition',pos4)



%A graficar por fin!!!
imp(s_est.tpo_predic:end) = imp(s_est.tpo_predic);
t_predic_a_fin = s_est.tpo_predic+1:length(I)-1;
soc_predic_a_fin = soc_counting(s_est.tpo_predic+1:end-1);


figure(2), sh = subplot(2,1,1); plot(soc_counting(1:s_est.tpo_predic),imp(1:s_est.tpo_predic),'k','Linewidth',2);%Aca deberia ir el real valor de R
%figure(2), sh = subplot(2,1,1), plot(R_soc_soc,R_soc(1:end-5),'b');
figure(2), subplot(2,1,1), hold on
figure(2), subplot(2,1,1), plot(soc_predic_a_fin,impedancia_media(2:end),'color','k','Linewidth',2,'Linestyle','--');
figure(2), subplot(2,1,1), plot([soc_counting(s_est.tpo_predic) soc_counting(s_est.tpo_predic)],[min(min(imp),min(lim_inf95_x1)) max(max(imp),max(lim_sup95_x1))],'k--')
figure(2), subplot(2,1,1), plot(soc_predic_a_fin,lim_inf95_x1(2:end),'color','k','Linewidth',1,'Linestyle','--');
figure(2), subplot(2,1,1), plot(soc_predic_a_fin,lim_sup95_x1(2:end),'color','k','Linewidth',1,'Linestyle','--');
figure(2), subplot(2,1,1), legend('Filtered x_1(|Z_{in}|)','Predicted x_1(|Z_{in}|)')
figure(2), subplot(2,1,1), axis tight


%figure(2), subplot(2,1,1), legend('Data |Z_{int}|','Filtered |Z_{int}|','Location','NorthWest');
figure(2), xlabel('SOC')
figure(2), ylabel('|Z_{int}| [\Omega]')
set(sh,'XDir','reverse');   
figure(2), sh = subplot(2,1,2); plot(soc_counting,smooth(V+I*s_est.CI(1),10),'k');
figure(2), subplot(2,1,2), hold on
figure(2), subplot(2,1,2), plot(soc_counting(1:s_est.tpo_predic),v_estim(1:s_est.tpo_predic)+I(1:s_est.tpo_predic)*s_est.CI(1),'color',naranjo,'Linewidth',2);
figure(2), subplot(2,1,2), plot(soc_counting(s_est.tpo_predic+1:end),s_result.reali(reali_plot).voltage+s_result.reali(reali_plot).Isint'*s_est.CI(1),'color',naranjo,'Linewidth',2,'Linestyle','--')
figure(2), subplot(2,1,2), axis([0 1 0 mod.VL+0.35])  
figure(2), subplot(2,1,2), plot([soc_counting(s_est.tpo_predic) soc_counting(s_est.tpo_predic)],[0 mod.VL+0.35],'k--')
figure(2), subplot(2,1,2), xlabel('SOC')
figure(2), ylabel('Voltage [V]')
figure(2), legend('V_{data}','V_{filtered}','V_{predicted}','Location','SouthWest')
set(sh,'XDir','reverse');
set(gca,'children',flipud(get(gca,'children')))



tof_real = find(soc_counting <0.05,1);
%display(s_result.mean.contNull)
%
h3 = figure(3);  figure(3), hold on;
figure(3), h3(1) = plot(soc_filtrado(1:s_est.tpo_predic)*100,'color',rojo,'Linewidth',2);
figure(3), h3(3) = plot(soc_counting(1:s_est.tpo_predic)*100,'k');
figure(3), h3(4) = plot(t_predic_a_fin,s_result.reali(reali_plot).soc(2:end)*100,'color','k','Linewidth',2,'Linestyle','--');

%figure(3), h3(5) = plot([prog.reali(reali_plot).max_IC prog.reali(reali_plot).max_IC]/mod.dt,[0 110],'color',verde_claro,'Linewidth',2,'Linestyle','--');
%figure(3), h3(6) = plot([prog.reali(reali_plot).min_IC prog.reali(reali_plot).min_IC]/mod.dt,[0 110],'color',verde_claro,'Linewidth',2,'Linestyle','--');
%figure(3), h3(7) = plot([prog.reali(reali_plot).tof prog.reali(reali_plot).tof],[0 110],'color',verde_osc,'Linewidth',2,'Linestyle','--');


%figure(3), h3(8) = plot(t_predic_a_fin,lim_inf(2:end)*100,'color',gris,'Linewidth',1,'Linestyle','--');
%figure(3), h3(9) = plot(t_predic_a_fin,lim_sup(2:end)*100,'color',gris,'Linewidth',1,'Linestyle','--');
figure(3), h3(10) = plot([s_est.tpo_predic s_est.tpo_predic],[0 110],'color','k','Linestyle','--');

figure(3), h3(11) = patch([0 length(soc_counting) length(soc_counting) 0],[4.5 4.5 5.5 5.5], morado_claro,'EdgeColor','none');
figure(3), h3(12) = plot(t_predic_a_fin,lim_inf95_x2(2:end)*100,'color','k','Linewidth',1,'Linestyle','--');
figure(3), h3(13) = plot(t_predic_a_fin,lim_sup95_x2(2:end)*100,'color','k','Linewidth',1,'Linestyle','--');
figure(3), h3(14) = stem(tof_real,5,'color','k','LineWidth',1);

%figure(3), set(h3(11), 'FaceAlpha', 0.65)
figure(3), ylim([0 105])
figure(3), xlim([0 length(I)*1.01])
figure(3), xlabel('Time [sec]')
figure(3), ylabel('SOC [%]')
figure(3), title('Particle Filters: Battery SOC')
figure(3), legend([h3(3) h3(1) h3(4) h3(11) h3(14)],'Offline SOC_{Ground Truth}','SOC_{filtered}','SOC_{predicted}','EOD zone','Ground Truth EOD','Location','NorthEast');



minIC = s_result.reali(reali_plot).min_IC;
maxIC = s_result.reali(reali_plot).max_IC;
tof = s_result.reali(reali_plot).tof;

figure(4), den = patch(t_predic_a_fin,s_result.reali(reali_plot).densidad(2:end)/max(s_result.reali(reali_plot).densidad),color_densidad);
figure(4), set(den, 'edgecolor',color_densidad-0.2)
%figure(4), stairs(t_predic_a_fin,prog.reali(reali_plot).densidad(2:end)/max(prog.reali(reali_plot).densidad),'color',color_densidad);
figure(4), hold on
figure(4), plot([maxIC maxIC]/mod.dt,[0 1.1],'color',verde_claro,'Linewidth',2,'Linestyle','--')
figure(4), plot([minIC minIC]/mod.dt,[0 1.1],'color',verde_claro,'Linewidth',2,'Linestyle','--')
figure(4), plot([tof tof],[0 1.1],'color',verde_osc,'Linewidth',2,'Linestyle','--')
figure(4), axis([minIC-3*(tof-minIC) maxIC+3*(maxIC-tof) 0 1]);
figure(4), xlabel('Time [sec]');
figure(4), title('Normalized SOC Density Function [sec]');
