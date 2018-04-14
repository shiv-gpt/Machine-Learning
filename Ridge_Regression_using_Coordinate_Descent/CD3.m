function [loss_val, lmda, W_val good_lmbda] = CD(feat, num_data, stopping_condition, lambda, X, Y, K)

    randn('state', 1);

    f = feat;
    n = num_data;
    s_c = stopping_condition;
    lambda_values = 10;
    
    %l = 2*sum(X*(Y - mean(Y)), 1);
    l = 2*max(X*(Y - mean(Y)));
    l1 = l
    %l
   
    A = 2*sum(X.*X, 2);
    C = zeros(f, 1);
    R = zeros(n, 1);
    W = zeros(f, 1);
    W_actual = zeros(f,1);
    W_actual(1:K, 1) = 10;
    %nnz(W_actual)
    B = 0;
    loss_delta = intmax('int64');
    loss = 0;
    precision = zeros(1, 10);
    recall = zeros(1,10);
    lmda = zeros(1,10);
    loss_val = zeros(1, 10);
    RMSE_t = zeros(1,10);
    RMSE_v = zeros(1, 10);
    NZ = zeros(1, 10);
    W_val = zeros(f, 10);
    count = 1;
    good_lmbda = -1*1000;
    while count <= 15%l >= 1
        count
        c = 1;
        loss_delta = intmax('int64');
        loss = 0;
        while (c < 300) && (loss_delta > 0.1)
            R = Y - X'*W - B;
            %R
            %pause(5);

            B_old = B;
            %B = mean(Y - X'*W);
            %B = mean(Y - X'*W - R);
            B = mean(R + B_old);
            R = R + B_old - B;
            for k = 1 : f
                %C(k,1) = 2*(sum( X(k,:)*R) + sum(X(k,:)*(X(k,:)'))*W(k,1));
                C(k,1) = 2*(sum( X(k,:)*R)) + A(k,1)*W(k,1);
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
            %if count == 10
                %disp('Here');
            %end
            %if(0 == mod(c,2))
                %C
                %pause(2);
            %end
            %if(0 == mod(c, 10))
         
            if count == 10
                disp('Here');
            end     %end
            loss_new = sum(R.*R) + l*sum(abs(W));
            loss_delta = abs(loss - loss_new);
            loss = loss_new;
            c = c + 1;
            loss_delta
            c
        end
        
        
        %W
        lmda(1, count) = l;
        [I, J, V] = find(W_actual);
        [I1, J1, V1] = find(W);
        %precision(1, count) = nnz(W(1:K,1))/nnz(W);
        %recall(1, count) = nnz(W(1:K,1))/K;
        loss_val(1, count) = loss;
        W_val(:, count) = W;
        
        RMSE_v(1, count) = Validate(W, B);
        if count > 1 && RMSE_v(1, count) > RMSE_v(1, count - 1)
            good_lmbda = lmda(1, count -1);
        end
        Output = X'*W + B;
        RMSE_t(1, count) = sqrt(sum((Y - Output).*(Y - Output))/size(W, 1));
        NZ(1, count) = nnz(W);
        count = count + 1;
        l = l/2;
        
        
    end
    W
    RMSE_v
    RMSE_t
    
    %precision
%     recall
    lmda
    Output = X'*W + B
    
    figure(1)
    plot(lmda, RMSE_t);
    title('RMSE - Training Plot');
    xlabel('Lambda');
    ylabel('RMSE_t');
    figure(2)
    plot(lmda, RMSE_v);
    title('RMSE - Validation Plot');
    xlabel('Lambda');
    ylabel('RMSE_v');
    figure(3)
    plot(lmda, NZ);
    title('Number of Non-Zeros Plot');
    xlabel('Lambda');
    ylabel('NNZ');
    
%     ylim([0 2]);
%     legend('Precision','Recall','Location','northeast');
end

