

%calculating alpha from francios et al (1982)
%check yourself at %check yourself at http://resource.npl.co.uk/acoustics/techguides/seaabsorption/
f = f0:f1; %the frequencies we care about
f = f/1000;%changing it to kHz
f = f';
D = 579; %depth of hydrophone in meters
S = 34.2552; %salinity at hydrophone depth from CalCofi
pH = 7.589; %pH at hydrophone depth from CalCofi
T = 6.5784; % temperature in celcius at depth of hydrophone from CalCOfi
c = 1412 + 3.21*T + 1.19*S +0.0167*D; % sound speed at hydrophone depth from calCofi



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

alpha = alpha';