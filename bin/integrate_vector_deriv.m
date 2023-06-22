function val = integrate_deriv( u, v, deriv_dimension, dx_vec, corner_vec, size_vec, m_vec, n_vec, mask )
  %PURPOSE:
  %This function transfers derivatives onto the wieght functions
  %
  %so we need a weight vector <w_y, -w_x> that is the curl of a scalar w. Then
  %\int dA <u_d,v_d> \cdot <w_y,  -w_x>  = \int dA (v_x - u_y) w  if w and its derivatives vanish on the boundary.
  %

  val =       integrate_deriv( v, 2, dx, corner_vec, size_vec, m_column, n_column, mask );
  val = val - integrate_deriv( u, 1, dx, corner_vec, size_vec, m_column, n_column, mask );

  m_vec2 = m_vec;
  n_vec2 = n_vec;
  m_vec2(deriv_dimension) = m_vec2(deriv_dimension) - 1;
  n_vec2(deriv_dimension) = n_vec2(deriv_dimension) - 1;

  val =     - m_vec(deriv_dimension)*integrate_data( data, dx_vec, corner_vec, size_vec, m_vec2, n_vec,  mask );
  val = val + n_vec(deriv_dimension)*integrate_data( data, dx_vec, corner_vec, size_vec, m_vec,  n_vec2, mask );

  mask_gradient = cell(1,3);
  [mask_gradient{1} mask_gradient{2} mask_gradient{3}] = gradient3d( mask, dx_vec );
  val = val - integrate_data( data, dx_vec, corner_vec, size_vec, m_vec, n_vec, mask_gradient{deriv_dimension} );
end
