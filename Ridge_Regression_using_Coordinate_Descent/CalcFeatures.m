function [Max_I, Min_I] = CalcFeatures(W)
    [Max, Max_I] = maxk(W, 10);
    [Min, Min_I] = mink(W, 10);
%     D = load('featureTypes.txt');
%     top = zeros(1, 10);
%     below = zeros(1, 10);
%     
%     for i =  1:10
%         top(1, i) = D(Max_I(1, i));
%         below(1, i) = D(Min_I(1, i));
%     end
end