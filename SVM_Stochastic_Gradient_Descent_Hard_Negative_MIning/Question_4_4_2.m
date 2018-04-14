function [alpha, HardNeg, prec, rec, Objective_function, valap] = Question_4_4_2(iterations, filepath)
    [trD, trLb, valD, valLb, trRegs, valRegs] = HW2_Utils.getPosAndRandomNeg();
    P_Ind = trLb == 1;
    PosD = trD(P_Ind);
    N_ind = trLb == -1;
    NegD = trD(N_ind);
    k = size(unique(trLb),1);
    num_epochs = 2000;
    C = 0.1;
    eta0 = 1;
    eta1 = 100;
%     [W, Accuracy_tr, Accuracy_val] = SVM_Q4(trD, trLb, valD, valLb, eta0, eta1, C, k, num_epochs);
%     [alpha, W, b, Result, Accuracy, trLb] = SVM_quadratic(10, 2, 0.1, trD, trLb);
    
    
    %Quadratic Programming
    [alpha, W, b] = SVM_quad_4_4_2(trD, trLb, C, 0);
%     size(alpha)
    
    threshold = 0.005;
%     W = trD*(alpha.*trLb);
%     I = find(((alpha<threshold)&(trLb==-1))|(trLb==1));
      I = (alpha > threshold & trLb == -1) | trLb == 1;
%     trD = trD(I);
    
    
    prev_best_ap = -inf;
    Objective_function = zeros(1, iterations);
    valap = zeros(1, iterations);
    for j = 1:iterations
%         trD = trD(I);
%         trLb = trLb(I);
        trD(:, I==0) = [];
        trLb(I==0, :) = [];
        HardNeg = HardestNeg(filepath, W, b, 0.5);
        Neg_Labels = -ones(size(HardNeg,2), 1);
        trD = cat(2, trD, HardNeg);
        trLb = cat(1, trLb, Neg_Labels);
%         size(trD)
%         size(trLb)
        [alpha, W, b] = SVM_quad_4_4_2(trD, trLb, C, 0);
        I = (alpha > threshold & trLb == -1) | trLb == 1;
        
        
        %%%%%
        Result = W'*trD + b;
        i = Result >= 1.0;
        Result(i) = 1;
        i = Result < 1.0;
        Result(i) = -1;

        Diff = trLb - Result';
        n = size(trLb, 1);
        Accuracy = ((n - nnz(Diff))/n)*100;
        Accuracy
        %%%%%
        
        
        %%%%%%%%%%
        Result = W'*valD + b;
        i = Result >= 1.0;
        Result(i) = 1;
        i = Result < 1.0;
        Result(i) = -1;

        Diff = valLb - Result';
        n = size(valLb, 1);
        Accuracy_val = ((n - nnz(Diff))/n)*100;
        Accuracy_val
        %%%%%%%%%%
        
        HW2_Utils.genRsltFile(W, b, 'val', '/home/shivang/Desktop/CSE_512/hw2_q4/hw2data/Output_4_4_3');
        [val_ap, prec, rec] = HW2_Utils.cmpAP('/home/shivang/Desktop/CSE_512/hw2_q4/hw2data/Output_4_4_3', 'val');
        j
        val_ap
        
        if(val_ap > prev_best_ap)
            prev_best_ap = val_ap;
            HW2_Utils.genRsltFile(W, b, 'test', '/home/shivang/Desktop/CSE_512/hw2_q4/hw2data/Output_4_4_3_test_best');
        end
        Objective_function(1,j) = sum(W.*W)/2 + C*sum(max(1 - valLb.*((W'*valD + b)'), 0));
        valap(1,j) =  val_ap;
    end
    
%     prec
%     rec(
    
    HW2_Utils.genRsltFile(W, b, 'test', '/home/shivang/Desktop/CSE_512/hw2_q4/hw2data/Output_4_4_3_test_final');
%     [test_ap, prec, rec] = HW2_Utils.cmpAP('/home/shivang/Desktop/CSE_512/hw2data/Output_4_4_3_test', 'test');
%     test_ap
%     prec
%     rec
    

end