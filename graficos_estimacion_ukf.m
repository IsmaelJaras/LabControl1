%Colores
rojo = [0.89, 0, 0.13];
morado_claro = [0.6, 0.4, 0.8];
violeta = [0.54, 0.17, 0.89];
naranjo  =  [0.91, 0.41, 0.17];
verde_claro  =  [0.13, 0.55, 0.13];
verde_osc =  [0, 0.27, 0.13];
gris = [0.3804    0.3804    0.3804];
color_densidad = [0.34 0.63 0.83];
azul_claro = [0, 0.5, 1];
rojo_claro = [0.88, 0.24, 0.19];


%Posicion de las figuras
fig2 = figure(2);
fig3 = figure(3);
scnsize = get(0,'ScreenSize');
pos2 = [0,scnsize(4) * (1/10),scnsize(3)/2,scnsize(4)*9/10];
pos3 = [scnsize(3)/2,pos2(2),pos2(3),pos2(4)];
set(fig2,'OuterPosition',pos2)
set(fig3,'OuterPosition',pos3)


figure(2), sh = subplot(2,1,1); plot(soc_counting(1:s_est.tpo_predic),imp(1:s_est.tpo_predic),'Color','k','Linewidth',2);%Aca deberia ir el real valor de R
figure(2), subplot(2,1,1), hold on
figure(2), subplot(2,1,1), plot(soc_counting(1:s_est.tpo_predic),imp(1:s_est.tpo_predic)+1.96*std_zin','--','color',azul_claro)
figure(2), subplot(2,1,1), plot(soc_counting(1:s_est.tpo_predic),imp(1:s_est.tpo_predic)-1.96*std_zin','--','color',azul_claro)
figure(2), subplot(2,1,1), legend('Filtered x_1(SOC)','Location','SouthEast')
figure(2), subplot(2,1,1), axis tight
ax_lim = axis;
figure(2), subplot(2,1,1), plot([soc_counting(s_est.tpo_predic) soc_counting(s_est.tpo_predic)],[ax_lim(3), ax_lim(4)],'k--')
xlim([0, 1])
figure(2), xlabel('SOC')
figure(2), ylabel('|Z_{int}| [\Omega]')
set(sh,'XDir','reverse');   
figure(2), sh = subplot(2,1,2); plot(soc_counting(1:s_est.tpo_predic),smooth(V(1:s_est.tpo_predic)+I(1:s_est.tpo_predic)*s_est.CI(1),10),'Color','k');
figure(2), subplot(2,1,2), hold on
figure(2), subplot(2,1,2), plot(soc_counting(1:s_est.tpo_predic),v_estim(1:s_est.tpo_predic)+I(1:s_est.tpo_predic)*s_est.CI(1),'Color',naranjo,'Linewidth',2,'Linestyle','--');
figure(2), subplot(2,1,2), axis([0 1 0 mod.VL+0.35])  
figure(2), subplot(2,1,2), plot([soc_counting(s_est.tpo_predic) soc_counting(s_est.tpo_predic)],[0 mod.VL+0.35],'k--')
figure(2), subplot(2,1,2), xlabel('SOC')
figure(2), ylabel('Voltage [V]')
figure(2), legend('V_{data}','V_{filtered}','Location','West')
set(sh,'XDir','reverse');
set(gca,'children',flipud(get(gca,'children')))


h3 = figure(3);  figure(3), hold on;
figure(3), plot(soc_filtrado(1:s_est.tpo_predic)*100,'Color','k','Linewidth',2);
figure(3), plot(soc_counting(1:s_est.tpo_predic)*100,'Color',violeta,'Linewidth',2);
figure(3), plot([s_est.tpo_predic s_est.tpo_predic],[0 110],'color','k','Linestyle','--');
figure(3), plot(soc_filtrado(1:s_est.tpo_predic)*100+196*std_soc,'--','color',rojo_claro);
figure(3), plot(soc_filtrado(1:s_est.tpo_predic)*100-196*std_soc,'--','color',rojo_claro);

figure(3), ylim([0 105])
figure(3), xlim([0 length(I)*1.01])
figure(3), xlabel('Time [sec]')
figure(3), ylabel('SOC [%]')
figure(3), title('Unscented Kalman Filter: Battery SOC')
figure(3), legend('SOC_{filtered}','Offline SOC_{Ground Truth}','Location','NorthEast');
set(gca,'children',flipud(get(gca,'children')))
