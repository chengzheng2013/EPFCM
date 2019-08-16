%  ---------------------------------------------------------------------------------------------------------------------------------------
%  General information
%  ---------------------------------------------------------------------------------------------------------------------------------------
%	This code provides IT2 EPFCM algorithm described in "Interval and General Type-2 Enhanced Possibilistic Fuzzy C-Means clustering"
%   which is appearing in Applied Soft Computing Journal.
%   Author: Shahabeddin Sotudian
%   Any kind of comments, suggestions, or bug reports are welcome and appreciated.
%   Please feel free to contact the author: sotudian AT bu DOT edu.
%% ---------------------------------------------------------------------------------------------------------------------------------------
function [V1,V2,UT1,UT2, E1,E2] = IT2_EPFCM_clustering (X, c,m1,m2,Theta1,Theta2,Cf,Cp,Eta)
Eta = Eta(:);
n = size(X, 1);
p = size(X, 2);
max_iter = 100;		 % Max. iteration
term_thr = 1e-4;     % Termination threshold
display = 1;		 % Display info or not
E1 = zeros(max_iter, 1);	% Array for termination measure values
E2 = zeros(max_iter, 1);	% Array for termination measure values
V1 = rand(c, p);
V2 = V1;
U1 = zeros (c, n);
T1 = zeros (c, n);
U2 = zeros (c, n);
T2 = zeros (c, n);

% Main loop
for i = 1:max_iter
 
    % fill the distance matrix
    dist1 = Distance_Function (V1, X);
    dist2 = Distance_Function (V2, X);
    % calculate new U
    tmp1 = dist1.^(-2/(m1-1));      
    U_new1 = tmp1./(ones(c, 1)*sum(tmp1));
    tmp2 = dist2.^(-2/(m2-1));      
    U_new2 = tmp2./(ones(c, 1)*sum(tmp2));
    U1=min(U_new1,U_new2);
    U2=max(U_new1,U_new2);

    % Correct the situation of "singularity" 
    si1 = find (tmp1 == Inf);
    U1(si1) = 1;
    si2 = find (tmp2 == Inf);
    U2(si2) = 1;
    if (size (si1, 1) ~= 0)
        disp('IT2 EPFCM, Warning: Singularity occured and corrected.');
    end
    if (size (si2, 1) ~= 0)
        disp('IT2 EPFCM, Warning: Singularity occured and corrected.');
    end

    % Claculate new T
    tmp1 = (Cp.*(dist1 .^ 2)) ./ ( Eta * ones (1, n));
    tmp1 = -tmp1;
    T_new1= nthroot(exp(tmp1),Theta1);
    tmp2 = (Cp.*(dist2 .^ 2)) ./ ( Eta * ones (1, n));
    tmp2 = -tmp2;
    T_new2= nthroot(exp(tmp2),Theta2);
    T1=min(T_new1,T_new2);
    T2=max(T_new1,T_new2);
    % Correct the situation of "singularity"
    T1(si1) = 1; 
    T2(si2) = 1;
    V_old1 = V1;
    V_old2 = V2;
    % Upper bound and Lower Bound Membership
    U_T1=(Cf.*U1)+(Cp.*T1);
    U_T2=(Cf.*U1)+(Cp.*T2);
    U_T3=(Cf.*U2)+(Cp.*T1);
    U_T4=(Cf.*U2)+(Cp.*T2);
    % Lower Bound Membership
    DUM1=min(U_T1,U_T2);
    DUM2=min(DUM1,U_T3);
    UT1=min(DUM2,U_T4);
    % Upper bound Membership
    DUM1=max(U_T1,U_T2);
    DUM2=max(DUM1,U_T3);
    UT2=max(DUM2,U_T4);
    % Karni-Mendel Alg
    for km=1:c
        UT_mean=(UT1+UT2)/2;
        maxUT_mean = max(UT_mean);

        index_c = find(UT_mean(km, :) == maxUT_mean);
        Y= X(index_c,:);

        a=UT1(km,index_c);
        b=UT2(km,index_c);
        a=a(:);
        b=b(:);
        F=cat(2,a,b);
        if size(Y,1)<=1
            a1=Y;
            a2=km;
            a3=index_c;
            a4=UT_mean;
            a5=maxUT_mean;
        end

        [XLeft,XRight,L,R]=KM_Alg(F,Y);

        V1(km,:)=XLeft';   
        V2(km,:)=XRight';
    end

    % check termination condition
    E1(i) = norm (V1 - V_old1, 1);
    E2(i) = norm (V2 - V_old2, 1);

        if display 
            fprintf('Iteration count = %d, Termination measure value => Lower= %f    Upper=%f\n', i, E1(i),E2(i));
        end

        if E1(i) <= term_thr, break; end,
        if E2(i) <= term_thr, break; end,
end

end