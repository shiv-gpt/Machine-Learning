function [alpha, W, b] = SVM_quad_4_4_2(trD, trLb, C, threshold)
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
end