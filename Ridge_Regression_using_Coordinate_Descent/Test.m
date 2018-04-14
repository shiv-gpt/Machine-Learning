function Ans =  Test(W, B)
    D = load('testData.txt');
    X = sparse(D(:,2), D(:,1), D(:,3));
    Output = X'*W + B;
    O = linspace(1, size(Output, 1), size(Output, 1));
    Ans = [O' Output];
    csvwrite('predTestLabels.csv',Ans);
end