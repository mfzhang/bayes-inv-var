clear all; clc;

nlpd = @(delta,sigma2) 0.5*(log(sigma2) +  (delta)^2*(1./sigma2));
sigma2 = 1e-4:0.001:100;


v_delta  = 0.1;
L = length(v_delta);
val = zeros(L,length(sigma2));
for i = 1 : L
    val(i,:) = nlpd(v_delta(i),sigma2);
end

semilogx(sigma2, val, 'LineWidth',2);
legend(mat2str(v_delta));
set(gca, 'FontSize', 14);
xlabel('\sigma^2');
ylabel('NLPD');
