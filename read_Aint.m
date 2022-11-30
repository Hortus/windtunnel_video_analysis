

bendfiles = dir('*bend*.txt');

height = csvread('windtunnel.csv');

oat_ind = [3,4,5,11,12,13,17,18,20];
wheat_ind = [1,2,6,7,8,9,10,14,15,16,19]; 

wheat_aveg_height = mean(height(wheat_ind,2))*10;% in mm 
oat_aveg_height = mean(height(oat_ind,2))*10; % in mm

for K = 1 : length(bendfiles)
  
  thisfilename = bendfiles(K).name;  %just the name
  bend_data = load(thisfilename); %load just this file
  area = bend_data(:,11);
  N=24;                           % number elements to average over
  mnX=mean(reshape(area(1:N*fix(numel(area)/N)),N,[]));
  Aint(K) = mnX(1,1);
end




wheat_aveg_Aint = mean(Aint(wheat_ind))*100; % in mm^2
oat_aveg_Aint = mean(Aint(oat_ind))*100; % in mm^2

rho_air = 1.225*10^(-6); %[g/mm^3]
%c_d = 0.61; %0.3; %drag air
c_d = 0.5; %requested

v_wind_oat = [6000:200:8400];
wind_load_oat = 1/2*rho_air*(oat_aveg_Aint)*c_d*v_wind_oat.^2;

wind_moment_oat = wind_load_oat*(wheat_aveg_height/2 - 47 -78);

