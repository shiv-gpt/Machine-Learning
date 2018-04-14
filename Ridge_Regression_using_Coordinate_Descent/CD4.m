function [RMSE_v, RMSE_t, Ans, W] = CD4(feat, num_data, lambda, X, Y)

    randn('state', 1);

    f = feat;
    n = num_data;
%     s_c = stopping_condition;
%     lambda_values = 10;
    
    %l = 2*sum(X*(Y - mean(Y)), 1);
    l = lambda;
%     l1 = l
    %l
   
    A = 2*sum(X.*X, 2);
    C = zeros(f, 1);
    R = zeros(n, 1);
    W = zeros(f, 1);
    top = 0;
    below = 0;
%     W_actual = zeros(f,1);
%     W_actual(1:K, 1) = 10;
    %nnz(W_actual)
    B = 0;
%     loss_delta = intmax('int64');
%     loss = 0;
%     precision = zeros(1, 10);
%     recall = zeros(1,10);
%     lmda = zeros(1,10);
%     loss_val = zeros(1, 10);
%     RMSE = zeros(1,10);
%     W_val = zeros(f, 10);
%     count = 1;
%     while count <= 15%l >= 1
%         count
        c = 1;
        loss_delta = intmax('int64');
        loss = 0;
        while (c < 300) && (loss_delta > 0.01)
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
         
%             if count == 10
%                 disp('Here');
%             end     %end
            loss_new = sum(R.*R) + l*sum(abs(W));
            loss_delta = abs(loss - loss_new);
            loss = loss_new;
            c = c + 1;
            loss_delta
            c
        end
        
        
        %W
%         lmda(1, count) = l;
%         [I, J, V] = find(W_actual);
%         [I1, J1, V1] = find(W);
        %precision(1, count) = nnz(W(1:K,1))/nnz(W);
        %recall(1, count) = nnz(W(1:K,1))/K;
        %loss_val(1, count) = loss;
        %W_val(:, count) = W;
        
        RMSE_v = Validate(W, B);
        Output = X'*W + B;
        RMSE_t = sqrt(sum((Y - Output).*(Y - Output))/size(W, 1));
        %[Top, Below] = CalcFeatures(W);
        Ans = Test(W, B);        
%         count = count + 1;
%         l = l/2;
        
        
    end
%     W
% %     figure
% %     precision
% %     recall
% %     lmda
%     Output = X'*W + B
%     plot(lmda, RMSE);
%     title('RMSE Plot');
%     xlabel('Lambda');
%     ylabel('RMSE');
%     ylim([0 2]);
%     legend('Precision','Recall','Location','northeast');


