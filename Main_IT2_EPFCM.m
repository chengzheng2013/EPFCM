%  ---------------------------------------------------------------------------------------------------------------------------------------
%  General information
%  ---------------------------------------------------------------------------------------------------------------------------------------
%	This code provides IT2 EPFCM algorithm described in "Interval and General Type-2 Enhanced Possibilistic Fuzzy C-Means clustering"
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
m1 = 2.0;
m2 = 4.0;
Theta1 = 3.0;
Theta2 = 5.0;
Cf=0.5;
Cp=0.5;
K=1;
% Initialization for IT2-EPFCM (Optional)
[V,U] = fcm(Xin,nC,[NaN 100 0.0001 0]);
ETA_init = Initialization_ETA (Xin, U, V, mean(m1,m2), K);
% IT2-EPFCM ------------------------------------------------------------
tic
[V1,V2,U1,U2, E1,E2] = IT2_EPFCM_clustering (Xin,nC,m1,m2,Theta1,Theta2,Cf,Cp,ETA_init);
toc
 
% Plotting -------------------------------------------------------------
% Plot input feature vectors
figure; plot(Xin(:,1),Xin(:,2),'o')
title ('Input feature vectors');
% Plot termination measure values
figure;
plot(E1);
hold on
plot(E2);
title ('Termination measure (IT2 EPFCM)');
xlabel ('Iteration num.');
ylabel ('Termination measure value');
cMarker = ['+' 'o' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h'];
cColor =  ['r' 'g' 'b' 'm' 'c' 'y' 'k' 'r' 'g' 'b' 'y' 'm' 'c'];
% Plot clustered feature vectors
figure;
maxU = max(U1);

cMarker = ['+' 'o' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h'];
cColor =  ['r' 'g' 'b' 'm' 'c' 'y' 'k' 'r' 'g' 'b' 'y' 'm' 'c'];

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

cMarker = ['+' 'o' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h'];
cColor =  ['r' 'g' 'b' 'm' 'c' 'y' 'k' 'r' 'g' 'b' 'y' 'm' 'c'];

for c = 1:nC
    index_c = find(U2(c, :) == maxU2);

    line(Xin(index_c, 1), Xin(index_c, 2), 'linestyle',...
        'none','marker', cMarker(c), 'color', cColor(c));
    
    hold on
    plot(V2(c,1),V2(c,2),['k' cMarker(c)],'markersize',15,'LineWidth',2)
end
title ('upper center');
% Plot membership functions
figure; hold on;
subplot (nC, 1, 1)
plot (U1(1, :), cColor(1))
hold on
plot (U2(1, :), cColor(1))
%title ('Membership functions');
for c = 2:nC
    subplot (nC, 1, c)
    plot (U1(c, :), cColor(c))
    hold on 
    plot (U2(c, :), cColor(c))
end

