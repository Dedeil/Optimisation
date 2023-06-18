%% 1.1 - Sous-differentielle et algorithme de sous-gradient
clear variables
close all
clc

sousGrad = @(x, alpha)(sign(x) .* (x~=0) + alpha .* (x==0.5));

lambda = .1;
xk = 3.27;
alpha = -1;
l = 50;
est = zeros(1, l);

for k=1:l
    tk = sousGrad(xk,alpha);
    xk = xk - lambda * tk;
    est(k) = xk;
end

figure
hold on
t = linspace(-3, 3, 3);
plot(t, abs(t));
plot(est, abs(est), '*-');

%% 1.2 - Enveloppe de Moreau
clear variables
close all
clc

a = 1; b = 0; gamma = .5;
f = @(x)(abs(a*x+b));
Mgamma = @(x)((a*x+b-a^2*gamma/2) .* (x>a*gamma-b/a) + ...
    (x-b/a).^2/(2*gamma) .* ((a*gamma-b/a>=x) & (x>=-a*gamma-b/a)) - ...
    (a*x+b+a^2*gamma/2 ).* (x<-a*gamma-b/a));

x = -5:10^-2:5;

figure
plot(x, Mgamma(x), x, f(x));

%% 1.3 - Operateur proximal
clear variables
close all
clc

gamma = .2; z = 1; x = -2:.1:2;
f = @(x)((x-z).^2);
y_chapeau = @(x)((2*gamma*z+x)/(2*gamma+1));

figure
plot(x, f(x), x, y_chapeau(x), x,x)
legend("f(x)=(x-z)²", "Operateur proximal")
grid on

%% 2.1 - Algorithme du point proximal
clear variables
close all
clc

gamma = 200;
epsilon = 10^-4;
argmin(1) = 5;
y_chapeau = @(x)((x-gamma) .* (x>gamma) + ...
    (x+gamma) .* (x<-gamma) + ...
    0 .* ((-gamma<=x) & (x<=gamma)));
dist = abs(argmin(1));
k = 1;
while epsilon < dist
    xk = argmin(k);
    argmin(k+1) = y_chapeau(xk);
    dist = sqrt(argmin(k)^2-argmin(k+1)^2);
    k = k + 1;
end

figure
t = linspace(-5, 5, 3);
hold on
plot(t, abs(t))
plot(argmin, abs(argmin), '*-');

%% 2.2 - Algorithme de gradient proximal (forward-backward)
clear variables
close all
clc

y_chapeau = @(x,gamma)((x-gamma) .* (x>gamma) + ...
    (x+gamma) .* (x<-gamma) + ...
    0 .* (-gamma<x<gamma)); 

lambda = 1;
f = @(x)(norm(x)^2);
gradf = @(x)(2*x);
g = @(x)(lambda*norm(x, 1));

gamma = .1;
epsilon = 10^-4;
n = 1;
argmin = randn(n, 1);
dist = norm(argmin(:, 1), 1);
fct_cout = [];
k = 1;

while epsilon < dist
    xk = argmin(:, k);
    argmin(:, k+1) = y_chapeau(xk-gamma*gradf(xk),lambda*gamma);
    dist = norm(xk-argmin(:, k+1), 2);
    fct_cout(k) = f(xk)+g(xk);
    k = k + 1;
end

figure
plot(1:k-1, fct_cout)

%% 3.1 - Résolution du modele du LASSO

%% 3.2 - Resolution du modele TV (Total Variation)
