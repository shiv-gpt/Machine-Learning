function [RMSE] = Validate(W, B)
    D = load('valData.txt');
    X = sparse(D(:,2), D(:,1), D(:,3));
    Y = load('valLabels.txt');
    Output = X'*W + B;
    RMSE = sqrt(sum((Y - Output).*(Y - Output))/size(W, 1));
end