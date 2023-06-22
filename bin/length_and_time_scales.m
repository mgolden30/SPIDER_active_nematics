function [length_scale, time_scale] = length_and_time_scales( Q11, Q12, mask, dx_vec )
  %{
  This calculates length and time scales based on the nematic tensor
  %}

  [Q11_x, Q11_y, Q11_t] = gradient3d(Q11, dx_vec);
  [Q12_x, Q12_y, Q12_t] = gradient3d(Q12, dx_vec);
  
  theta_time_deriv  = sqrt( Q11_t.^2 + Q12_t.^2 );
  theta_space_deriv = sqrt( Q11_x.^2 + Q12_x.^2 + Q11_y.^2 + Q12_y.^2 );
  
  %Don't trust numerical derivatives when the mask is almost zero.
  acceptable = mask > 0.1;
  
  length_scale = pi * sum(acceptable, 'all') / sum( theta_space_deriv.*acceptable, 'all' );
  time_scale   = pi * sum(acceptable, 'all') / sum( theta_time_deriv.*acceptable, 'all' );
end