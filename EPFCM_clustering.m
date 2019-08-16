%  ---------------------------------------------------------------------------------------------------------------------------------------
%  General information
%  ---------------------------------------------------------------------------------------------------------------------------------------
%	This code provides EPFCM algorithm described in "Interval and General Type-2 Enhanced Possibilistic Fuzzy C-Means clustering"
%   which is appearing in Applied Soft Computing Journal.
%   Author: Shahabeddin Sotudian
%   Any kind of comments, suggestions, or bug reports are welcome and appreciated.
%   Please feel free to contact the author: sotudian AT bu DOT edu.
%% ---------------------------------------------------------------------------------------------------------------------------------------


function [V, U, T, E,ObjFun_EPFCM] = EPFCM_clustering (X, c,m,Theta,Cf,Cp,Eta)
Eta = Eta(:);
n = size(X, 1);
p = size(X, 2);
max_iter = 100;		    % Max. iteration
term_thr = 1e-4;	    % Termination threshold
display = 1;		    % Display info or not
E = zeros(max_iter, 1);	% Array for termination measure values
ObjFun_EPFCM=zeros(max_iter, 1);
V = rand(c, p);
U = zeros(c, n);
T = zeros(c, n);

% Main loop  
for i = 1:max_iter,

    % fill the distance matrix
    dist = Distance_Function (V, X); 
    % calculate new U, suppose m != 1
    tmp = dist.^(-2/(m-1));      
    U = tmp./(ones(c, 1)*sum(tmp));
    % Correct the situation of "singularity".
    si = find (tmp == Inf);
    U(si) = 1;
    if (size (si, 1) ~= 0)
        display ('EPFCM, Warning: Singularity occured and corrected.');
    end
    % Claculate new T
    tmp = (Cp.*(dist .^ 2)) ./ ( Eta * ones (1, n));
    tmp = -tmp;
    T = nthroot(exp(tmp),Theta);
    T(si) = 1; % Correct the situation of "singularity".

    %   objective function 
    PU1 = U.^m;
    PU2 = T.^Theta;
    ObjFun_EPFCM(i) = Cf*(sum(sum((dist.^2).*PU1)))  +  Cp*sum(sum((dist.^2).*PU2))  +  sum(Eta.*sum((PU2.*log(PU2)-PU2),2));
    
    % new centers
    V_old = V;
    Us = Cf.*(U.^m);
    Ts = Cp.*(T.^Theta);
    V = ((Us+Ts)*X) ./ ((ones(p, 1)*sum((Us+Ts)'))'); 

    
    % check termination condition
    E(i) = norm (V - V_old, 1);
        if display, 
            fprintf('Iteration count = %d, Termination measure value = %f\n', i, E(i));
        end
    
        if E(i) <= term_thr
            break;
        end
end
