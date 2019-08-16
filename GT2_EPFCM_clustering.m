%  ---------------------------------------------------------------------------------------------------------------------------------------
%  General information
%  ---------------------------------------------------------------------------------------------------------------------------------------
%	This code provides GT2 EPFCM algorithm described in "Interval and General Type-2 Enhanced Possibilistic Fuzzy C-Means clustering"
%   which is appearing in Applied Soft Computing Journal.
%   Author: Shahabeddin Sotudian
%   Any kind of comments, suggestions, or bug reports are welcome and appreciated.
%   Please feel free to contact the author: sotudian AT bu DOT edu.
%% ---------------------------------------------------------------------------------------------------------------------------------------
function [V1,V2,UT1,UT2, E1,E2] = GT2_EPFCM_clustering (X,c,Eta)
Eta = Eta(:);
Cf=0.5;
Cp=0.5;
NAP=5;                      % Number of alfa plane

n = size(X, 1);
p = size(X, 2);
max_iter = 100;		        % Max. iteration
term_thr = 1e-4;		    % Termination threshold
display = 1;		        % Display info or not
E1 = zeros(max_iter, 1);	% Array for termination measure values
E2 = zeros(max_iter, 1);	% Array for termination measure values

V1 = rand(c, p);
V2 = V1;

U1 = zeros (c, n);
T1 = zeros (c, n);
U2 = zeros (c, n);
T2 = zeros (c, n);

%  Fuzzifier Gaussian  ----------------------------------------------
% Gaussian  Number for m
sigma1=3;
mean1=10;
% Gaussian  Number for t
sigma2=4;
mean2=15;
alfa=[];
counter=1;
for i=(1/NAP):(1/NAP):1
    alfa(counter)=i;
    % alfa cut
   M1(counter)=(mean1-sigma1*sqrt((-2)*log(alfa(counter))));
   M2(counter)=(mean1+sigma1*sqrt((-2)*log(alfa(counter))));
   Theta1(counter)=(mean2-sigma2*sqrt((-2)*log(alfa(counter))));
   Theta2(counter)=(mean2+sigma2*sqrt((-2)*log(alfa(counter))));
counter=1+counter;
end
% Plot Fuzzifiers
l=-5:0.5:(max(mean1,mean2)+4*max(sigma1,sigma2));
subplot(1,2,1)
h1=gaussmf(l,[sigma1 mean1]);
plot(l,h1,'c','linewidth',1.5)
hold on 
scatter(cat(2,M1,M2),cat(2,alfa,alfa),'b','LineWidth',4)
title('Fuzzifier for (m)')
axis([-5,(max(mean1,mean2)+5*max(sigma1,sigma2)),0,1])

subplot(1,2,2)
h2=gaussmf(l,[sigma2 mean2]);
plot(l,h2,'r','linewidth',1.5)
hold on
scatter(cat(2,Theta1,Theta2),cat(2,alfa,alfa),'b','LineWidth',4)
title('Fuzzifier for (t)')
axis([-5,(max(mean1,mean2)+5*max(sigma1,sigma2)),0,1])

 
% Main loop
for i = 1:max_iter,
    
    V_old1 = V1;
    V_old2 = V2;
     for con=1:NAP   
            [V1,V2, UT1_alfacut,UT2_alfacut ] = OneStep_GT2EPFCM(V1,V2,X,c,M1(con),M2(con),Theta1(con),Theta2(con),Cf,Cp,Eta ) ;
            CENTERS1(:,:,con)=V1;
            CENTERS2(:,:,con)=V2;
            U_matrixes1(:,:,con)=UT1_alfacut;
            U_matrixes2(:,:,con)=UT2_alfacut;
     end
     UT1=U_matrixes1;
     UT2=U_matrixes2;
     
    % Union Center
     T2_Center_Lower=zeros(size(CENTERS2(:,:,1)));
     T2_Center_Upper=zeros(size(CENTERS2(:,:,1)));
    for con=1:NAP
    T2_Center_Lower=T2_Center_Lower+(CENTERS1(:,:,con)*alfa(con));
    T2_Center_Upper=T2_Center_Upper+(CENTERS2(:,:,con)*alfa(con));
    end

    V1=T2_Center_Lower/sum(alfa);
    V2=T2_Center_Upper/sum(alfa);


    E1(i) = norm (V1 - V_old1, 1);
    E2(i) = norm (V2 - V_old2, 1);

    if display, 
        fprintf('Iteration count = %d, Termination measure value => Lower= %f    Upper=%f\n', i, E1(i),E2(i));
    end
    % check termination condition
    if E1(i) <= term_thr, break; end,
    if E2(i) <= term_thr, break; end,
end




end