%  ---------------------------------------------------------------------------------------------------------------------------------------
%  General information
%  ---------------------------------------------------------------------------------------------------------------------------------------
%	This code provides GT2 EPFCM algorithm described in "Interval and General Type-2 Enhanced Possibilistic Fuzzy C-Means clustering"
%   which is appearing in Applied Soft Computing Journal.
%   Author: Shahabeddin Sotudian
%   Any kind of comments, suggestions, or bug reports are welcome and appreciated.
%   Please feel free to contact the author: sotudian AT bu DOT edu.
%% ---------------------------------------------------------------------------------------------------------------------------------------

clc;
clear all;
close all;

% Data -----------------------------------------------------------------
load fcmdata.dat
Xin = fcmdata;
% Options --------------------------------------------------------------
nC = 2 ;  % Number of clusters
K=1;

% Initialization for GT2-EPFCM (Optional)
[V,U] = fcm(Xin,nC,[NaN 100 0.0001 0]);
m1=2;
m2=5;
Eta_int = Initialization_ETA (Xin, U, V, mean(m1,m2), K);

% GT2-EPFCM ------------------------------------------------------------
tic
[V1,V2,U_matrixes1,U_matrixes2, E1,E2] = GT2_EPFCM_clustering (Xin,nC,Eta_int);
toc
  
% Plotting -------------------------------------------------------------  
% Plot membership functions
cMarker = ['+' 'o' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h'];
cColor =  ['r' 'g' 'b' 'm' 'c' 'y' 'k' 'r' 'g' 'b' 'y' 'm' 'c'];
figure; 

for c = 1:nC
   subplot (nC, 1, c) 
    hold on;
plot (U_matrixes1(c,:,1),'r')
hold on
plot (U_matrixes2(c,:,1),'r')
hold on
plot (U_matrixes1(c,:,2),'b')
hold on
plot (U_matrixes2(c,:,2),'b')
hold on
plot (U_matrixes1(c,:,3),'g')
hold on
plot (U_matrixes2(c,:,3),'g')
hold on
plot (U_matrixes1(c,:,4),'m')
hold on
plot (U_matrixes2(c,:,4),'m')
    
end
 
% Plot termination measure values
figure;
plot(E1);
hold on
plot(E2);
title ('Termination measure (GT2 EPFCM)');
xlabel ('Iteration num.');
ylabel ('Termination measure value');

% Plot clustered feature vectors
figure;
U1=sum(U_matrixes1,3)/size(U_matrixes1,3);
U2=sum(U_matrixes2,3)/size(U_matrixes2,3);
  
maxU = max(U1);

for c = 1:nC
    index_c = find(U1(c, :) == maxU);

    line(Xin(index_c, 1), Xin(index_c, 2), 'linestyle',...
        'none','marker', cMarker(c), 'color', cColor(c));
    
    hold on
    plot(V1(c,1),V1(c,2),['k' cMarker(c)],'markersize',15,'LineWidth',2)
end
title ('lower center');

figure;
maxU2 = max(U2);

for c = 1:nC
    index_c = find(U2(c, :) == maxU2);

    line(Xin(index_c, 1), Xin(index_c, 2), 'linestyle',...
        'none','marker', cMarker(c), 'color', cColor(c));
    
    hold on
    plot(V2(c,1),V2(c,2),['k' cMarker(c)],'markersize',15,'LineWidth',2)
end
title ('upper center');

