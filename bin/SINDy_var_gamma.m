function [Xi, eta] = SINDy_var_gamma (Theta, dXdt, gamma)
%{
Compute sparse regression on dX = Theta * Xi
Regression technique used: sequential least squares

Modified procedure from:
S. H. Rudy, S. L. Brunton, J. L. Proctor, J. N. Kutz, Data-driven 
discovery of partial differential equations. Sci. Adv. 3, e1602614 (2017)
%}

%Xi = Theta \ dXdt;
Xi = pinv(Theta)*dXdt;
lambda = gamma*mean(abs(dXdt)); % threshold to determine term as "small"
eta = 0;

for i = 1:20
  product = zeros(size(Xi)); 
  [~,w] = size(Theta);
  for p_ind = 1:w
    product(p_ind) = Xi(p_ind)*mean(abs(Theta(:,p_ind)));
  end

  smallinds = (abs(product) < lambda);
  Xi(smallinds) = 0;                        % set negligible terms to 0
  for ind = 1:size(dXdt,2)   
    biginds = ~smallinds(:,ind);
    Xi(biginds,ind) = Theta(:,biginds) \ dXdt(:,ind);
  end

   %Switch to Patrick's notation
   Q   = Theta;
   q_0 = dXdt;
   c   = Xi;

   eta = norm(Q*c - q_0)/max(vecnorm(c'.*Q, 2, 1));
   fprintf('Iteration %d: eta = %f\n', i, eta);
end %i

end %function definition
