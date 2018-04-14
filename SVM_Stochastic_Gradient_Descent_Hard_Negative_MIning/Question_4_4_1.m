function [ap, prec, rec] =Question_4_4_1()
    [trD, trLb, valD, valLb, trRegs, valRegs] = HW2_Utils.getPosAndRandomNeg();
    k = size(unique(trLb),1);
    num_epochs = 2000;
    C = 10;
    eta0 = 1;
    eta1 = 100;
%     [W, Accuracy_tr, Accuracy_val] = SVM_Q4(trD, trLb, valD, valLb, eta0, eta1, C, k, num_epochs);
%     b = 0.1;
%     Accuracy_tr
%     Accuracy_val
    [alpha, W, b] = SVM_quad_4_4_2(trD, trLb, C, 0);
    HW2_Utils.genRsltFile(W, b, 'val', '/home/shivang/Desktop/CSE_512/hw2data/Output_4_4_1');
    [ap, prec, rec] = HW2_Utils.cmpAP('/home/shivang/Desktop/CSE_512/hw2data/Output_4_4_1', 'val');
    ap
    prec
    rec
end