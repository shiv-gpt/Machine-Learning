function [W, Accuracy_tr, Accuracy_val] = SVM_Q4(trD, trLb, valD, valLb, eta0, eta1, C, k, num_epochs)
%     load('/home/shivang/Desktop/CSE_512/hw2data/q3_1_data.mat');
    %load('/home/shivang/Desktop/CSE_512/hw2data/q3_2_data.mat');
%     size(trD)
%     size(trLb)
    [d, n] = size(trD);
    W = zeros(d, k);
    L = zeros(1,n);
    delta_W_yi = zeros(d,1);
    delta_W_y_i_hat = zeros(d,1);
    y_i = 1;
    y_i_hat = 1;
    
    if k == 2 
        t = trLb == -1;
        trLb(t) = 2;
    end
    
    for e = 1:num_epochs
        e
        eta = eta0/(eta1 + e);
        p = randperm(n);
        for i = 1:n
%             i
%             p(1,i)
            y_i = trLb(p(1,i), 1);
            [m,ind] = maxk(W'*trD(:,p(1,i)), 2);
            if ind(1,1) == y_i
                y_i_hat = ind(2,1);
            else
                y_i_hat = ind(1,1);
            end
%             y_i
%             y_i_hat+
            Loss = max(W(:,y_i_hat)'*trD(:,p(1,i)) - W(:,y_i)'*trD(:,p(1,i)) + 1, 0);
            del_Loss = C*trD(:,p(1,i));
            for j = 1:k
%                 W(:,1)'*trD(:,2)
%                 W(:,y_i)'*trD(:,p(1,i))
%                 Loss = max(W(:,y_i_hat)'*trD(:,p(1,i)) - W(:,y_i)'*trD(:,p(1,i)) + 1, 0);
%                 Loss
                if j == y_i
%                     size(delta_W_yi)
                    if Loss == 0
                        delta_W_yi = W(:,y_i)/n;
                    else
%                         delta_W_yi = W(:,y_i)/n - C*trD(:,p(1,i));
                          delta_W_yi = W(:,y_i)/n - del_Loss;
                    end
%                     size(trD(:,p(1,i)))
%                     size(delta_W_yi)
                    %W_old = W(:,y_i);
                    W(:,y_i) =W(:,y_i) - eta*delta_W_yi;
                elseif j == y_i_hat
                    if Loss == 0
                        delta_W_y_i_hat = W(:,y_i_hat)/n;
                    else
%                         delta_W_y_i_hat = W(:,y_i_hat)/n + C*trD(:,p(1,i));
                        delta_W_y_i_hat = W(:,y_i_hat)/n + C*trD(:,p(1,i));
                    end
%                     size(delta_W_y_i_hat)
%                      W_old = W(:,y_i_hat);
                    W(:,y_i_hat) = W(:,y_i_hat) - eta*delta_W_y_i_hat;
                else
                     W(:,j) = W(:,j) - eta*(W(:,j)/n);
                end                   
                
                        
                    
            end
        end
            
    end
    
    Output = W'*trD;
%     size(Output)
    [m, ind] = maxk(Output,1);
%     size(ind)
    Test = ind - trLb';
%     nnz(Test)
    Output_vald = W'*valD;
    [m, ind] = maxk(Output_vald,1);
    if k == 2
     t = valLb == -1;
     valLb(t) = 2;
    end
    Test_val = ind - valLb';
    Accuracy_tr = ((n - nnz(Test))/n)*100;
    Accuracy_val = ((n - nnz(Test_val))/n)*100;
end