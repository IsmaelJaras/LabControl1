# LabControl1
Archivos y contenido:

1."datos_bat1_robot.m": carga base de datos de robot móvil y parámetros de modelo asociados.

2."datos_E1_FUDS.m": carga base de datos FUDS y parámetros de modelo asociados.

3."E1D2_FUDS.mat": base de datos FUDS.

4."E1D3_identificacion2.mat": base de datos para identificación de la impedancia interna de la batería (utilizado por datos_E1_FUDS.m).

5."estimacion2v2.m": módulo de estimación basado en FP con modelo nuevo.

6."estimacion2v2_ofcl.m": módulo de estimación basado en FP con modelo nuevo y OFCL.

7."estimacion_ukf.m": módulo de estimación basado en UKF.

8."estimacion_ukf_ofcl5.m": módulo de estimación basado en UKF con OFCL.

9."graficos_bn.m": script para generar gráficos en B&N de los esquemas de estimación y pronóstico (sin partículas o intervalo de confianza en estimación).

10."graficos_color.m": script para generar gráficos a color de los esquemas de estimación y pronóstico (sin partículas o intervalo de confianza en estimación).

11."graficos_color_pf.m": script para generar gráficos a color del esquema de estimación con FP y pronóstico (gráfico de valor individual de cada partícula).

12."graficos_color_pf_ofcl.m": script para generar gráficos a color del esquema de estimación con FP y OFCL y pronóstico (gráfico de valor individual de cada partícula y activación del OFCL).

13."graficos_color_ukf.m": script para generar gráficos a color del esquema de estimación con UKF y pronóstico (gráfico de intervalo de confianza en estimación).

14."graficos_color_ukf_ofcl.m": script para generar gráficos a color del esquema de estimación con UKF y pronóstico (gráfico de intervalo de confianza en estimación y activación del OFCL).

15."graficos_estimacion.m": script para generar gráficos a color sólo de los módulos de estimación (sin partículas o intervalo de confianza en estimación).

16."graficos_estimacion_pf.m": script para generar gráficos a color sólo del módulo de estimación con FP (gráfico de valor individual de cada partícula).

17."graficos_estimacion_pf_ofcl.m": script para generar gráficos a color sólo del módulo de estimación con FP y OFCL (gráfico de valor individual de cada partícula y activación del OFCL).

18."graficos_estimacion_ukf.m": script para generar gráficos a color sólo del módulo de estimación con UKF (gráfico de intervalo de confianza en estimación).

19."graficos_estimacion_ukf_ofcl.m": script para generar gráficos a color sólo del módulo de estimación con UKF (gráfico de intervalo de confianza en estimación y activación del OFCL).

20."Intervalo.m": función para determinar el valor superior e inferior del intervalo de 95% de confianza del ToF pronosticado.

21."main_ukf.m": código principal para ejecución del esquema de estimación y pronóstico con estimación basada en UKF (con o sin OFCL).

22."main_v2.m": código principal para ejecución del esquema de estimación y pronóstico con estimación basada en FP (con o sin OFCL).

23."markov_I.m": función para generar valores de corriente futura en base a una cadena de Markov de dos estados.

24."PDF.m": función para calcular la pdf a partir de una cdf.

25."prediccion.m": función para realización individual de la predicción del SOC en un horizonte temporal dado.

26."pronostico.m": función principal del módulo de ponóstico en base FP.

27."prueba_bat1_Irobot.mat": base de datos robot móvil.

28."regularizacion.m": función de regularización para FP en pronóstico.

29."resample_sys_dp.m": función de resampling para FP.
