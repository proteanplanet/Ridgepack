
cont=[28:2:62];

cnames={'EC_60_30','EC_60_30_Degraded','EC_60_30_Old'};

for i=1:length(cnames)

coastname=char(cnames{i}); % grid names

centlat=-40; % degrees north
centlon=180; % degrees east
horizon=10; % degrees of satellite horizon (0-90)
ridgepack_E3SM_resolution_check(centlat,centlon,horizon,coastname,cont)

centlat=0; % degrees north
centlon=180; % degrees east
horizon=5; % degrees of satellite horizon (0-90)
ridgepack_E3SM_resolution_check(centlat,centlon,horizon,coastname,cont)

centlat=40; % degrees north
centlon=180; % degrees east
horizon=10; % degrees of satellite horizon (0-90)
ridgepack_E3SM_resolution_check(centlat,centlon,horizon,coastname,cont)

centlat=-40; % degrees north
centlon=-25; % degrees east
horizon=10; % degrees of satellite horizon (0-90)
ridgepack_E3SM_resolution_check(centlat,centlon,horizon,coastname,cont)

centlat=0; % degrees north
centlon=-25; % degrees east
horizon=5; % degrees of satellite horizon (0-90)
ridgepack_E3SM_resolution_check(centlat,centlon,horizon,coastname,cont)

centlat=28; % degrees north
centlon=-72; % degrees east
horizon=5; % degrees of satellite horizon (0-90)
ridgepack_E3SM_resolution_check(centlat,centlon,horizon,coastname,cont)

centlat=90; % degrees north
centlon=-45; % degrees east
horizon=5; % degrees of satellite horizon (0-90)
ridgepack_E3SM_resolution_check(centlat,centlon,horizon,coastname,cont)

centlat=58; % degrees north
centlon=-50; % degrees east
horizon=5; % degrees of satellite horizon (0-90)
ridgepack_E3SM_resolution_check(centlat,centlon,horizon,coastname,cont)

centlat=25; % degrees north
centlon=-80; % degrees east
horizon=5; % degrees of satellite horizon (0-90)
ridgepack_E3SM_resolution_check(centlat,centlon,horizon,coastname,cont)

centlat=67; % degrees north
centlon=-168; % degrees east
horizon=5; % degrees of satellite horizon (0-90)
ridgepack_E3SM_resolution_check(centlat,centlon,horizon,coastname,cont)

end






