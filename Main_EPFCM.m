%  ---------------------------------------------------------------------------------------------------------------------------------------
%  General information
%  ---------------------------------------------------------------------------------------------------------------------------------------
%	This code provides EPFCM algorithm described in "Interval and General Type-2 Enhanced Possibilistic Fuzzy C-Means clustering"
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
m = 2;
Theta = 3;
Cf=0.5;
Cp=0.5;
K=1;
% Initialization for EPFCM (Optional)
[V,U] = fcm(Xin,nC,[NaN 100 0.0001 0]);
ETA_init = Initialization_ETA (Xin, U, V, m, K);

% EPFCM ----------------------------------------------------------------
tic
[V_EPFCM,U_EPFCM,T_EPFCM,E,ObjFun_EPFCM] =EPFCM_clustering (Xin,nC,m,Theta,Cf,Cp,ETA_init);
toc
% Plotting -------------------------------------------------------------
% Plot input feature vectors
figure; plot(Xin(:,1),Xin(:,2),'o')  
title ('Input feature vectors');
% Plot termination measure values
figure;
plot(E);
title ('Termination measure (EPFCM)');
xlabel ('Iteration num.');
ylabel ('Termination measure value');
% Plot clustered feature vectors
cMarker = ['+' 'o' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h'];
cColor =  ['r' 'g' 'b' 'm' 'c' 'y' 'k' 'r' 'g' 'b' 'y' 'm' 'c'];
figure;
maxU = max(U);
for c = 1:nC
    index_c = find(U(c, :) == maxU);

    line(Xin(index_c, 1), Xin(index_c, 2), 'linestyle',...
        'none','marker', cMarker(c), 'color', cColor(c));
    
    hold on
    plot(V(c,1),V(c,2),['k' cMarker(c)],'markersize',15,'LineWidth',2)
end
title ('Clustered feature vectors (EPFCM)');
% Plot membership functions
figure; hold on;
subplot (nC, 1, 1)
plot (U(1, :), cColor(1))
title ('Membership functions (EPFCM)');
for c = 2:nC
    subplot (nC, 1, c)
    plot (U(c, :), cColor(c))
end


