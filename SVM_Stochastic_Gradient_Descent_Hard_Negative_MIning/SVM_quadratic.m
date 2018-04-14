function [alpha, W, b, Result, Accuracy, trLb, Objective_function, confusion, Support_Vectors] = SVM_quadratic(C, k, threshold, trD, trLb)
%     load(filepath);
    
    H  = (trLb*trLb').*(trD'*trD);
    H = cast(H, 'double');
    f = -ones(1, size(trLb, 1));
    A = zeros(1, size(trLb, 1));
    b = 0;
    Aeq = trLb';
    Aeq = cast(Aeq, 'double');
    beq = 0;
    lb = zeros(size(trLb, 1), 1);
    ub = C*ones(size(trLb, 1),1);
    alpha = quadprog(H,f,A,b,Aeq,beq,lb,ub);
    W = trD*(alpha.*trLb);
    I = min(find((alpha>threshold)&(trLb==1)));
    g = trD'*trD;
    b = 1 - g(I,:)*(alpha.*trLb);
    
    Result = W'*trD + b;
    i = Result >= 1.0;
    Result(i) = 1;
    i = Result < 1.0;
    Result(i) = -1;
    
    Diff = trLb - Result';
    n = size(trLb, 1);
    Accuracy = ((n - nnz(Diff))/n)*100;
    
    
    
    
    load('/home/shivang/Desktop/CSE_512/hw2data/q3_1_data.mat');
    W
    Result = W'*valD + b;
    Result
    i = Result >= 1.0;
    Result(i) = 1;
    i = Result < 1.0;
    Result(i) = -1;
    
    Diff = valLb - Result';
    n = size(valLb, 1);
    Val_Accuracy = ((n - nnz(Diff))/n)*100;
    Val_Accuracy
   
    Objective_function = sum(W.*W)/2 + C*sum(max(1 - valLb.*((W'*valD + b)'), 0));
%     Objective_function = sum(W.*W)/2 + C*sum(max(1 - valLb.*((Result)'), 0));
     confusion = confusionmat(valLb, Result');
     
    Objective_function
    NS = alpha > 0.001*C;
    Support_Vectors = alpha(NS);
    size(Support_Vectors);
%     P = trD*trD';
   
    
%      [d, n] = trD;
%     
%     H = zeros(d+1, d+1
end