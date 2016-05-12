function AirDensity = cjw_airdensity(Pressure,Temperature)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute air density from pressure and temperature
%simple relationship only, for extra precision use an external climatology!
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%10/JAN/2014
%
%inputs
%---------
%
%Temperature - temperature  (K)
%Pressure - pressure (hPa)
%
%outputs
%---------
%
%AirDensity - air density, kg/m^3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Rd = 287.; %Gas constant for dry air

AirDensity = Pressure .* 1e2  ./ (Rd .* Temperature);


