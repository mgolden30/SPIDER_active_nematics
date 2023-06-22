function [c_ave, residual_ave, residual_sequence_ave, singular_ave, c_std, residual_std, residual_sequence_std, singular_std, models] = svd_regression_statistics( A, threshold, subsamples, size_samples )
  %{
   This does svd regression on an ensemble of samples to do statistics
  %}

  num_samples = subsamples; %Number of times you perform SVD regression for statistics
  
  [nw, nl] = size(A);
  if( nw < size_samples)
    fprintf('Error: number of samples is larger than the number of integration domains.\n')
    return;
  end
  
  models = 0;
  cs                 = [];
  residuals          = []; %Just a scalar
  residual_sequences = [];
  singular_vals      = [];
  for sample=1:num_samples
      y = randsample(nw, size_samples);
      small_A = A(y,:);      
      [c, residual, residual_sequence] = svd_regression(small_A, threshold);
      
      if( c(1) == 1 )
        models = models + 1;
      end
      %if( c(1) == 1 )
        cs = [cs c];
        residuals = [residuals residual];
        residual_sequences = [ residual_sequences residual_sequence];
      %end
  end
  
  %Take the transpose of everything
  cs
  scatter(1:subsamples, vecnorm(A*cs), 'filled')
  
  cs = cs';
  residuals = residuals';
  residual_sequences = residual_sequences';
  singular_vals = singular_vals';
  
  
  %Now just take mean and std
  c_ave = mean(cs);
  c_std =  std(cs);
  
  residual_ave = mean(residuals);
  residual_std =  std(residuals);
  
  residual_sequence_ave = mean(residual_sequences);
  residual_sequence_std =  std(residual_sequences);
  
  singular_ave = mean(singular_vals);
  singular_std =  std(singular_vals);
end