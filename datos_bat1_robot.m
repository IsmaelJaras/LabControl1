load prueba_bat1_Irobot;

soc_real     = 1-cumsum(V.*I)./sum(V.*I);

l_datos = 2920;
I = I(1+data_offset:l_datos);
V = V(1+data_offset:l_datos);



r0 = 0.3;

%====================================================
%                Parámetros del Modelo
%====================================================

mod.dt = 1;
mod.Ecrit = 20934;

%ruidos del mod
%mod.std_x = [0.0011 0.004526];
mod.std_x = [0.0011 0.0045];
mod.std_obs = 0.057;
%500


mod.VL = 3.9773;
mod.beta = 16;
mod.alfa = 0.0776;
mod.g = 19.6445;
%mod.g2 = 7.0088;
mod.Vo = 4.12;
