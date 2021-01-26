function [TL_SS] = range2SS_TL(f0, f1, r, zs, zr, pH, T, S, c)

% take Range, Depth, Frequency, Salinity, Temperature, and pH to
% compute the frequency dependent spherical spreading between a source and
% a receiver
%
% f0 = start frequency (Hz)
% f1 = end frequency (Hz)
% zr = depth of receiver (m)
% zs = depth of source (m)
% r = horizontal range to receiver (m)
% pH = pH
% T = temperature
% S = salinity
% c = sound speed (m/s)

f = f0:f1; %range of frequencies
f = f/1000; %changing it to kHz
f = f';
D = zr

slantRange = sqrt((r^2) + ((zs-zr)^2))

%Boric Acid Contribution
A1 = (8.86/c)*10^((0.78*pH) -5);
P1 = 1;
F1 = (2.8*((S/35)^0.5))*10^(4-(1245/(273+T)));

%MgSO4 contribution
A2 = 21.44*(S/c)*(1+0.025*T);
P2 = 1-((1.37*(10^-4))*D) + (6.2*(10^-9))*D^2;
F2 = (8.17*10^(8-1990/(273+T)))/(1+0.0018*(S-35));

%Pure Water Contribution
A3 = (4.937*10^-4) - (2.59*10^-5)*T + (9.11*(10^-7))*T^2 - (1.50*(10^-8))*T^3;
P3 = 1-(3.83*(10^-5))*D + (4.9*(10^-10))*D^2;


%alpha equation
%total absorption = boric acid contrib. + MgSO4 contrib. + P. Water contrib.
for n = 1:length(f)
    alpha(n) = (A1*P1*F1*(f(n)^2))/(f(n)^2 + F1^2) + (A2*P2*F2*(f(n)^2))/((f(n)^2)+F2^2) +A3*P3*(f(n)^2);
end
%%%the units for alpha are dB/km, need to get it in meters like our slant
%%%range
alpha = alpha/1000;
alpha = alpha';
slantRange = slantRange'
TL = 20*log10(slantRange);%spherical spreading transmission loss at each frequency

%adding alpha to the spherical spreading model for frequency dependent TL
TL_SS = TL + alpha*(slantRange);
end

