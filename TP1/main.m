%% Partie 1
clear variables
close all

lambda = .4;
epsilon = 10^-4;
xk = [5, 12];
df = @(x)(2 * x);

d = length(xk);
if d==1
        figure
        hold on
        t = linspace(xk+1, -(xk+1), 1000);
        plot(t, t.^2)
end

[xk, grad, i] = descenteGradient(xk, lambda, epsilon, df);

%% Partie 2.1
clear variables
close all

addpath("matrices", "operators")

t = linspace(0, 2 * pi, 100);
H = matH(size(t),'gaussian',3);
xbar = sin(t)';
n = randn(size(xbar)) * .05;
z = H * xbar + n;

% Modele des moindres carres
xMoindresCarres = (H' * H)\(H' * z);

% Methode de descente de gradient approche
lambda = .4;
epsilon = 10 ^ -5;
df = @(x)(2 * H' * (H * x - z));
[xDescenteGradient] = descenteGradient(zeros(100, 1), lambda, epsilon, df);

% Modele de Tikhonov 1D approche
lambdaT = .7;
op = 'i';
s = size(t);
switch op
    case 'i'
        Gamma = eye(s(2));
    case 'g'
        Gamma = matGamma(s, 'gradient');
    case 'l'
        Gamma = matGamma(s, 'laplacian');
end

df = @(x)(2 * H' * (H * x - z) + 2 * lambdaT * (Gamma' * Gamma) * x);
[xTik1] = descenteGradient(z, lambda, epsilon, df);

% Modele de Tikhonov 1D exact
xTikExact = (H' * H + lambdaT * (Gamma' * Gamma)) \ (H' * z);

figure
hold on
plot(t, xbar)
plot(t, z, "g")
plot(t, xMoindresCarres, "--")
plot(t, xDescenteGradient)
plot(t, xTik1, 'r')
plot(t, xTikExact, "b--")

legend("xbar", "z", "xMoindresCarres", "xDescenteGradient", "xTik", "xTikExact")
xlim([-inf, inf])

%% Partie 2.3
clear variables
close all

addpath("matrices", "operators")

img = im2double(imread("cameraman.tif"));

n = randn(size(img)) * .2;
z = opH(img, "gaussian", 5) + n;

figure
subplot(231)
imshow(img)
subplot(232)
imshow(z)

% Modele de Tikhonov 2D approche
lambdaT = .7;
lambda = .4;
epsilon = 10 ^ -5;
imgTik = zeros(size(z));
dImg = @(x)(2 * opHt(opH(x, "gaussian", 5) - z, "gaussian", 5) + 2 * lambdaT * opDt(opD(x)));
fcost = @(x)(norm(opH(x, "gaussian", 5)) ^ 2 + lambdaT * norm(opD(x)) ^ 2);
grad = dImg(imgTik);

%TowDcost = fcost(z);

i = 0;
while norm(grad) > epsilon
    imgTik = imgTik - lambda * grad;
    grad = dImg(imgTik);
    i = i + 1;
    %TowDcost = [TowDcost, fcost(imgTik)];
end

subplot(233)
imshow(imgTik)
subplot(2, 3, [4, 5, 6])
%plot(TowDcost)

mse = norm(imgTik - img, 'fro') ^ 2;
snr = 10 * log10(norm(img, 2) ^ 2 / mse);
