function test(W)
    load('/home/shivang/Desktop/CSE_512/hw2data/q3_2_data.mat');
    %load('/home/shivang/Desktop/CSE_512/HW2/weights.mat');
    Output = W'*tstD;
    [m, ind] = maxk(Output,1);
%     size(ind)
    O = linspace(1, size(ind, 2), size(ind, 2));
    Ans = [O' ind'];
    csvwrite('predTestLabels.csv',Ans);
end