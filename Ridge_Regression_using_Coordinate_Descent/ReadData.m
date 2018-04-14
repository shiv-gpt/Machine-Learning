function [X, Y] = ReadData()
    D = load('trainData.txt');
    X = sparse(D(:,2), D(:,1), D(:,3));
    Y = load('trainLabels.txt');
end