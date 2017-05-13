Falias = Broadband_Tools.getAliasingFrequency(SYS.Main_Setup)/2/pi*SYS.signal_info.c / 1e3;
Falias = min(Falias,1.5938);
SimCircBright=Res_Matrix{1}(Hrz_Vec>0.15 & Hrz_Vec<Falias);
SimCircQuiet=Res_Matrix{2}(Hrz_Vec>0.15 & Hrz_Vec<Falias);
RealCircQuiet=Res_Matrix{4}(Hrz_Vec>0.15 & Hrz_Vec<Falias);
RealCircBright=Res_Matrix{3}(Hrz_Vec>0.15 & Hrz_Vec<Falias);
mean(RealCircQuiet)
mean(SimCircQuiet)
mean(RealCircBright*2-RealCircQuiet*2)
mean(SimCircBright*2-SimCircQuiet*2)