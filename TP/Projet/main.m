%% Chargement du model
clear variables
close all
clc

addpath("fonctions", "images/3D-denoising", "matrices")

[xbar, faces] = loadOff("teapot.off");
%% Affichage
close all
figure
subplot(221)
trisurf(faces, xbar(:,1), xbar(:,3), xbar(:,2))
set(gca,'DataAspectRatio',[1,1,1])

%% Bruitage du model
subplot(222)
z = xbar + .05 * randn(size(xbar));
trisurf(faces, z(:,1), z(:,3), z(:,2))
set(gca,'DataAspectRatio',[1,1,1])

%% Algorithme de debruitage
lambdaT = .6;
lambda = .01;
epsilon = 10^-1;
xchapeau = z;
% grad = @(x)(2 * (x - z + lambdaT * opDt(opD(x))));
% R = grad(xchapeau);
% i = 0;
% fcost = [];
% while norm(R) > epsilon
%     xchapeau = xchapeau - lambda * R;
%     R = grad(xchapeau);
%     i = i + 1;
%     fcost = [fcost, norm(xchapeau)];
% end

op = 'l';
switch op
    case 'g'
        R = @(x)(matGamma3D(x, faces, 'gradient'));
    case 'l'
        R = @(x)(matGamma3D(x, faces, 'laplacian'));
end
fcost = [];
df = @(x)(2 * (x - z + lambdaT * R(x) * x));
grad = df(xchapeau);
i = 0;
while norm(grad) > epsilon
    xchapeau = xchapeau - lambda * grad;
    grad = df(xchapeau);
    fcost = [fcost, norm(xchapeau)];
    i = i + 1;
%     if i > 50
%         break
%     end
end

subplot(223)
trisurf(faces, xchapeau(:,1), xchapeau(:,3), xchapeau(:,2))
set(gca,'DataAspectRatio',[1,1,1])

subplot(224)
plot(fcost)
