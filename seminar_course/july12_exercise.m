filename = "TGDTGair.txt";
air = importdata(filename);
filename = "TGDTGn2.txt";
n2 = importdata(filename);

Tair = air.data(:,5);
Tn2 = n2.data(:,5);
DTGair = air.data(:,2);
DTGn2 = n2.data(:,2);
TGair = air.data(:,3);
TGn2 = n2.data(:,3);

%air
figure
plot(Tair,TGair,'r')
title("TG VS T")
figure
plot(Tair,DTGair,'b')
title("DTG VS T")

%n2
figure
plot(Tn2,TGn2,'r')
title("TG VS T")
figure
plot(Tn2,DTGn2,'b')
title("DTG VS T")