
%Colores
naranjo  =  [1.0000    0.7000         0];
verde_claro  =  [0.4    0.4         0.4];
verde_osc =  [0.2000    0.2000         0.2];
morado_claro =  [0.6    0.6    0.6];
gris = [0.3804    0.3804    0.3804];
color_densidad = [0.85 0.85 0.85];



%Posicion de las figuras
fig2 = figure(2);
fig3 = figure(3);
scnsize = get(0,'ScreenSize');
borders = get(fig2,'OuterPosition') - get(fig2,'Position');
edge = -borders(1)/2;
pos2 = [edge,scnsize(4) * (1/4),scnsize(3)/2 - edge,scnsize(4)*3/4];
pos3 = [scnsize(3)/2 + edge,pos2(2),pos2(3)*0.92,pos2(4)/2];
set(fig2,'OuterPosition',pos2)
set(fig3,'OuterPosition',pos3)


figure(2), sh = subplot(2,1,1); plot(soc_counting(1:s_est.tpo_predic),imp(1:s_est.tpo_predic),'Color','b','Linewidth',2);%Aca deberia ir el real valor de R
figure(2), subplot(2,1,1), hold on
figure(2), subplot(2,1,1), plot([0,0],[min(imp),min(imp)])

figure(2), subplot(2,1,1), plot([soc_counting(s_est.tpo_predic) soc_counting(s_est.tpo_predic)],[min(imp), max(imp)],'k--')
figure(2), subplot(2,1,1), legend('Filtered x_1(SOC)')
figure(2), subplot(2,1,1), axis tight


figure(2), xlabel('SOC')
figure(2), ylabel('|Z_{int}| [\Omega]')
set(sh,'XDir','reverse');   
figure(2), sh = subplot(2,1,2); plot(soc_counting(1:s_est.tpo_predic),smooth(V(1:s_est.tpo_predic)+I(1:s_est.tpo_predic)*s_est.CI(1),10),'Color','k');
figure(2), subplot(2,1,2), hold on
figure(2), subplot(2,1,2), plot(soc_counting(1:s_est.tpo_predic),v_estim(1:s_est.tpo_predic)+I(1:s_est.tpo_predic)*s_est.CI(1),'Color','g','Linewidth',2,'Linestyle','--');
figure(2), subplot(2,1,2), axis([0 1 0 mod.VL+0.35])  
figure(2), subplot(2,1,2), plot([soc_counting(s_est.tpo_predic) soc_counting(s_est.tpo_predic)],[0 mod.VL+0.35],'k--')
figure(2), subplot(2,1,2), xlabel('SOC')
figure(2), ylabel('Voltage [V]')
figure(2), legend('V_{data}','V_{filtered}','Location','NorthEast')
set(sh,'XDir','reverse');


h3 = figure(3);  figure(3), hold on;
figure(3), plot(soc_filtrado(1:s_est.tpo_predic)*100,'Color','r','Linewidth',2);
figure(3), plot(soc_counting(1:s_est.tpo_predic)*100,'k');
figure(3), plot([s_est.tpo_predic s_est.tpo_predic],[0 110],'color','k','Linestyle','--');

figure(3), ylim([0 105])
figure(3), xlim([0 length(I)*1.01])
figure(3), xlabel('Time [sec]')
figure(3), ylabel('SOC [%]')
figure(3), title('Particle Filters: Battery SOC')
figure(3), legend('SOC_{filtered}','Offline SOC_{Ground Truth}','Location','NorthEast');
