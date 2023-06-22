function [c_ave, res_ave, res_seq_ave, c_std, res_std, res_seq_std] = ...
         model_discovery( A, threshold, times_to_subsample, size_of_subsamples )
  %{
   PURPOSE:
   The purpose of this function is to take a matrix of integrated library
   terms "A" and find as many models as possible.

   INPUT:
   A         - a matrix of integrated library terms
   threshold - see documentation of svd_regression
   alignment - determines how close two values are. alignment ∈ [0,1].
   error     - determines how close a value is to zero:
        Suppose we have two models:
            model_1: "1.0000 0 0 0 0   "
            model_2: "1.0000 0 0 0 0.01"
        if you think the last elements of model_1 and model_2 are too
        close, and thus model_1 and model_2 should be the same model, you
        just set "error" larger than "0.01" (last element of model_2). Then
        the regression should average model_1 and model_2:
            model_new: "1.0000 0 0 0 0.005"


   OUTPUT:
   G_new - appends "A" and "c_ave" so that we can run regression on "G_new"
           to see if we could find more new models.
   c_ave - a matrix of average coefficients.
  %}
  
  [nw, nl] = size(A);
  assert( size_of_subsamples < nw );
  
  cs                 = zeros( nl, times_to_subsample );
  residuals          = zeros( 1,  times_to_subsample );
  residual_sequences = zeros( nl, times_to_subsample );
  
  for sample=1:times_to_subsample
      y = randsample(nw, size_of_subsamples); small_A = A(y,:);
      [c, residual, residual_sequence] = svd_regression( small_A, threshold );
      %c 
      if sample > 1
        c1 = cs(:,1);
        if c'*c1 < 0
          c = -c; 
        end
      end
      
      %Make sure everything is a column vector.
      assert( size(c,1)         ~= 1 );
      assert( size(residual_sequence,1) ~= 1 );
      
      cs(:,sample) = c;
      residuals(sample) = residual;
      residual_sequences(:,sample) = residual_sequence;
  end
  
  %This plots the residual for each model discovered.
  %scatter(1:times_to_subsample, vecnorm(A*cs), 'filled')
  
  my_ones = 0*residuals + 1; %Need for std function.
  
  c_ave = mean( cs, 2 );
  c_std = std( cs, my_ones, 2 );
  res_ave = mean(residuals);
  res_seq_ave = mean(residual_sequences, 2); 
  res_std = std(residuals);
  res_seq_std = std(residual_sequences, my_ones,2);
  
  %So I can see what model discovery picks out
  cs
end