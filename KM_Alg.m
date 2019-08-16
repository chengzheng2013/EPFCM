function [XLeft,XRight,L,R]=KM_Alg(F,X)

%% KM Algorithm for Computing X Left

% a) Sort X matrix and Membership
Norm_X = X;
for i=1:size(Norm_X,1)
   lowerX(i,:)= norm(Norm_X(i,:));
end
    
[~, ind] = sort(lowerX);
Sorted_X = Norm_X(ind,:);
sorted_F = F(ind,:);

% b) Initializeation
Teta=(sorted_F(:,1)+sorted_F(:,2))/2;
Sorted_X=Sorted_X';
isZero=(sum(Teta)==0);
if isZero
    C_prime=0;
else
    C_prime=(Sorted_X*Teta)/sum(Teta);
end

Counter=0;
while(1)
    % c) Find switch point k (1 <= k <= N ? 1) such that Xk <= X <= Xk+1
    sPointLeft = 0;
    for i=1:size(Norm_X,1)-1;
        if norm(C_prime)>=norm(Sorted_X(:,i)) && norm(C_prime)<=norm(Sorted_X(:,i+1))
            sPointLeft = i;
            break
        end
    end
%    a1= sPointLeft
        
    % d) Compute C_left
    for i=1:size(Norm_X,1)
        if i<=sPointLeft
            fn(i,1) = sorted_F(i,2);
        elseif i>sPointLeft
            fn(i,1) = sorted_F(i,1);
        end
    end
    
    if(sum(fn)==0)
        C_left=0;
    else
        C_left = (Sorted_X*fn)/sum(fn);
    end
%     a2=C_left
%     a3=abs(norm(C_prime)-norm(C_left))
   % e) if C_prime==C_left stop else go to c)
    if(abs(norm(C_prime)-norm(C_left))<10^-3)
        
        XLeft = C_prime;
        L = sPointLeft;
        break;
    else
        Counter=Counter+1;
        if Counter>=(size(Norm_X,1)/2)
            if(abs(norm(C_prime)-norm(C_left))<10^-2)
                XLeft = C_prime;
                L = sPointLeft;
                break;
            end
        end
        C_prime=C_left;
    end
end


%% KM Algorithm for Computing X Right

% a) Sort X matrix and Membership
Norm_X = X;
for i=1:size(Norm_X,1)
   lowerX(i,:)= norm(Norm_X(i,:));
end
    
[~, ind] = sort(lowerX);
Sorted_X = Norm_X(ind,:);
sorted_F = F(ind,:);

% b) Initializeation
Teta=(sorted_F(:,1)+sorted_F(:,2))/2;
Sorted_X=Sorted_X';
isZero=(sum(Teta)==0);
if isZero
    C_prime=0;
else
    C_prime=(Sorted_X*Teta)/sum(Teta);
end

Counter=0;
while(1)
    % c) Find switch point k (1 <= k <= N ? 1) such that Xk <= X <= Xk+1
    sPointRight = 0;
    for i=1:size(Norm_X,1)-1;
        if norm(C_prime)>=norm(Sorted_X(:,i)) && norm(C_prime)<=norm(Sorted_X(:,i+1))
            sPointRight = i;
            break
        end
    end
    
        
    % d) Compute C_Right
    for i=1:size(Norm_X,1)
        if i<=sPointRight
            fn_right(i,1) = sorted_F(i,1);
        elseif i>sPointRight
            fn_right(i,1) = sorted_F(i,2);
        end
    end
    
    if(sum(fn_right)==0)
        C_right=0;
    else
        C_right = (Sorted_X*fn_right)/sum(fn_right);
    end
    
    
    
    
   % e) if C_prime==C_right stop else go to c)
    if(abs(norm(C_prime)-norm(C_right))<10^-3)
        
        XRight = C_prime;
        R = sPointRight;
        break;
    else
        Counter=Counter+1;
        if Counter>=(size(Norm_X,1)/2)
            if(abs(norm(C_prime)-norm(C_right))<10^-2)
                XRight = C_prime;
                R = sPointRight;
                break;
            end
        end
        C_prime=C_right;
    end
end

