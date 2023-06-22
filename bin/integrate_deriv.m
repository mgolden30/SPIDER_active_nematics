function val = integrate_deriv( data, deriv_dimension, dx_vec, corner_vec, size_vec, m_vec, n_vec, mask )
  %PURPOSE:
  %This function transfers derivatives onto the wieght functions
  %val = \int dV \partial_i data * weight = - \int dV data \partial_i weight  

  dd = deriv_dimension;
  m2 = m_vec;
  n2 = n_vec;
  m2(dd) = m2(dd) - 1;
  n2(dd) = n2(dd) - 1;

  val =     - m_vec(dd)*integrate_data( data, dx_vec, corner_vec, size_vec, m2, n_vec,  mask );
  val = val + n_vec(dd)*integrate_data( data, dx_vec, corner_vec, size_vec, m_vec,  n2, mask );

  persistent mask_gradient;
  if isempty(mask_gradient)
    mask_gradient = cell(1,3);
    [mask_gradient{1}, mask_gradient{2}, mask_gradient{3}] = gradient3d( mask, dx_vec );
  end

  val = val - integrate_data( data, dx_vec, corner_vec, size_vec, m_vec, n_vec, mask_gradient{deriv_dimension} );
end
