function mask = make_mask( theta, dx, alpha, beta )
  %{
  PURPOSE:
  Finds the topological defects of the nematic field. This is done by looking at the determinant of the gradient of
  the director field. This was suggested by Roman
 
  IN:
  Q_11 - component of the nematic tensor
  Q_12 - "
  dx - grid spacing
  alpha and beta - parameters to finely tune the mask.
                   I designed them so that alpha=beta=1 should be a good
                   default choice.

  boolean_mask_branch_cuts - a boolean (true/false) that determines if the
  branch cuts where n changes sign should be masked.

  OUT:
  mask - a 2+1D matrix
  %}
  
  n1 = cos(theta);
  n2 = sin(theta);

  %Calculate gradient
  dx_vec = [dx dx dx];
  if boolean_mask_branch_cuts
    %If this is true, we want to kill branch cut regions. Take the usual
    %finite differences.
    [nx_y, nx_x, ~] = gradient3d( nx, dx_vec );
    [ny_y, ny_x, ~] = gradient3d( ny, dx_vec );
  else
    %Check alignment of neighbors so that branch cuts pose no difficulty.
    [nx_x, nx_y, ~, ny_x, ny_y, ~] = nematic_derivative(nx, ny, dx_vec);
  end
  
  det = nx_x.*ny_y - nx_y.*ny_x;
  
  nan_check( det, 'det grad n' );
  
  width = 2;
  mask = imgaussfilt3( abs(det), width ); %Blur this so that it is smooth
  
  mean_log = sum( log(mask), 'all' )/gpx/gpy/gpt;
  nan_check( mean_log, 'mean_log' );
  
  %use 1e-9 to make sure there is no division by zero.
  mask = tanh( alpha*(exp(mean_log))./(mask + 1e-9) );
  mask = mask.^beta;
  
  mask = imgaussfilt3( mask, width/2 );
  
  nan_check( mask, 'mask in the make_mask function' );
end