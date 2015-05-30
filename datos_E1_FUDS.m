load E1D3_identificacion2.mat;
Vi = V;
Ii = I;
soc_iden = 1-cumsum(Vi.*Ii)./sum(Vi.*Ii);

load E1D2_FUDS.mat;

soc_real = 1-cumsum(V.*I)./sum(V.*I);

l_datos = length(V);
I = I(1+data_offset:l_datos);
V = V(1+data_offset:l_datos);



%r0 = 0.10;%era 0.12 pero se aproximo
r0 = 0.10;

%====================================================
%                Parámetros del Modelo
%====================================================

mod.dt = 1;
mod.Ecrit = 46858;


%ruidos del mod
%mod.std_x = [0.0015 0.0055];
mod.std_x = [0.0015 0.0055];
mod.std_obs = 0.067;
%500


mod.VL = 3.9974;
mod.beta = 17;
mod.alfa = 0.1469;
mod.g = 10.4954;
%mod.g2 = 7.0088;
mod.Vo = 4.144;
