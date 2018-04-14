function [l1, C, W, A] = CD(feat, num_data, stopping_condition, lambda, X, Y, K)

    randn('state', 1);

    f = feat;
    n = num_data;
    s_c = stopping_condition;
    
    %l = 2*sum(X*(Y - mean(Y)), 1);
    l = 2*max(X*(Y - mean(Y)));
    l1 = l
    l
   
    A = 2*sum(X.*X, 2);
    C = zeros(f, 1);
    R = zeros(n, 1);
    W = zeros(f, 1);
    %W(1:K, 1) = 10;
    B = 0;
    for c = 1:300
        R = Y - X'*W - B;
        %R
        %pause(5);
        
        B_old = B;
        %B = mean(Y - X'*W);
        B = mean(Y - X'*W - R);
        R = R + B_old - B;
        for k = 1 : f
            C(k,1) = 2*(sum( X(k,:)*R) + sum(X(k,:)*(X(k,:)'))*W(k,1));
            W_old = W;
            if(C(k,1) < -1*l)
                W(k, 1) = (C(k,1) + l)/A(k,1);
            elseif(C(k,1) >= -1*l && C(k,1) <= l)
                W(k,1) = 0;
            else
                W(k,1) = (C(k,1) - l)/A(k,1);
            end
            R = R + X'*W_old - X'*W;        
                
               
        end
        %if(0 == mod(c,2))
            %C
            %pause(2);
        %end
        if(0 == mod(c, 10))
            l = l/2;
        end
       
        
    end
    W
    
end

