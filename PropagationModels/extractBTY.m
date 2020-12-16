function [R, bath] = extractBTY(slat, slon, hlat, hlon)
% take ship lat/lon and hydrophone lat/lon and make bathymetry file from
% between the two
% NOTE: I'm using the hydrophone as the source and ship as receiver for now


%% load in bathymetry data
load('sbc_bathymetry')

%% interpolate w/ interp2 to extract a line of bathymetry?
lati = linspace(min([slat, hlat]), max([slat, hlat]), 100);
loni = linspace(min([slon, hlon]), max([slon, hlon]), 100);

[xi, yi] = latlon2xy(lati, loni, hlat, hlon);

ri = sqrt(xi.^2 + yi.^2);

zi = interp2(lat, lon, z, lati, loni);

bath = abs(zi);
R = ri;
%% Plot (for testing)
% If you want to verify that this extracted the correct section of the
% bathymetry, uncomment all the plotting, and put in the hydrophone depth
% and source depth here:

% hz = -565;
% sz = -10
% 
% Nlat = find(lat>hlat-.09 & lat<slat+.09);
% Nlon = find(lon>hlon-.09 & lon<slon+.09);
% 
% 
% 
% plot3(loni, lati, bathi, 'linewidth', 1.5)
% hold on
% plot3([hlon hlon], [hlat hlat], [-abs(hz) -abs(sz)], 'r^:')
% plot3([slon slon], [slat slat], [-abs(hz) -abs(sz)], 'rx:')
% surf( lon(Nlon), lat(Nlat), z(Nlon, Nlat).')
% colormap jet
% shading interp
% hold off
% title('Extracted bathymetry between source and receiver')
% legend('interpolated bathy line', 'hydrophone', 'ship location')