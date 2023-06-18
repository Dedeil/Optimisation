function [x, grad, i] = descenteGradient(x, lambda, epsilon, df)
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
    grad = df(x);
    i = 0;
    while norm(grad) > epsilon
        x = x - lambda * grad;
        grad = df(x);
        i = i + 1;
    end
end
