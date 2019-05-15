function [alim,blim,rhot] = coefcdf(mat,mbt,Caat,Cbbt,Cabt)
alim=-mat./sqrt(Caat); blim=-mbt./sqrt(Cbbt);
rhot = Cabt./sqrt(Caat.*Cbbt);
end