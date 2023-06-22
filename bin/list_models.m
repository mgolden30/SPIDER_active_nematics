function [c, model_string] = list_models(G, labels, threshold)
  %{
  PURPOSE:
  The purpose of this function is to perform regression and output a LaTeX
  string for ease of use.

  INPUT:
  G - integrated matrix.
  labels - a cell full of latex strings for each term.

  OUTPUT:
  c - coefficient vector
  model_string - Latex string
  %}
  
  [~, nl] = size(G); %nl = num_library

  [c, residual, residual_sequence] = svd_regression( G, threshold );
  
  m  = (c ~= 0);
  num_m = sum(m);

  labels_sparse = cell(num_m, 1);
  j=1;
  for i=1:numel(c)
    if c(i) ~= 0 
      labels_sparse{j} = labels{i};
      j = j+1;
    end
  end
  
  %If single term model is found, just output it.
  if num_m == 1
    model_string = sprintf("NA") + " & $" + labels_sparse{1} + "$\\";
    return;
  end
  
  %Check for uncertainty in the model
  G2 = G(:,m); %restrict ourselves to the proper subspace
  [c_ave, res_ave, res_seq_ave, c_std, res_std, res_seq_std] = model_discovery( G2, threshold, 500, size(G,[1])/2 );

  outstring =  sprintf("$%.3e \\pm %.1f", res_ave, res_std/res_ave*100) + "\% $ & $";
  for i = 1:sum(m)
    outstring = outstring + sprintf(" + (%.3f \\pm %.1f", c_ave(i), c_std(i)/abs(c_ave(i))*100 ) + "\% )" + labels_sparse{i}; 
  end
  outstring = outstring + "$\\";
  
  model_string = outstring;
  c = zeros(nl,1);
  c(m) = c_ave;
end