function val = integrate_second_deriv( data, dd1, dd2, dx_vec, corner_vec, size_vec, m_vec, n_vec, mask )
  %{
  PURPOSE:
  This function transfers 2 derivatives of the same coordinate onto the wieght functions
  val = \int dV \partial_i \partial_j data * weight = \int dV data \partial^2_i weight  
  
  IN:
  data - spacetime data you want to integrate
  dd1  - derivative dimension 1
  dd2  - derivative dimension 2
  %}
  
  m = m_vec; %shorthand
  n = n_vec; %shorthand
  
  %save derivatives of the mask so we don;t have to recompute them
  %This is a big speedup, but it means you need to clear memeory when
  %switching data sets.
  persistent mask_gradient;
  if isempty(mask_gradient)
    mask_gradient = cell(1,3);
    [mask_gradient{1} mask_gradient{2} mask_gradient{3}] = gradient3d( mask, dx_vec );
  end

  persistent mask_double_gradient;
  if isempty(mask_double_gradient)
    mask_double_gradient = cell(3,3);
    [mask_double_gradient{1,1} mask_double_gradient{1,2} mask_double_gradient{1,3}] = gradient3d( mask_gradient{1}, dx_vec );
    [mask_double_gradient{2,1} mask_double_gradient{2,2} mask_double_gradient{2,3}] = gradient3d( mask_gradient{2}, dx_vec );
    [mask_double_gradient{3,1} mask_double_gradient{3,2} mask_double_gradient{3,3}] = gradient3d( mask_gradient{3}, dx_vec );
  end
  
  %All derivatives on mask
  val = integrate_data( data, dx_vec, corner_vec, size_vec, m, n, mask_double_gradient{dd1,dd2} );
  
  
  % Now all derivatives on polynomial
  m2 = m;      %copy m
  a  = m2(dd1); %save what we have to reduce
  m2(dd1) = m2(dd1)-1;
  b = m2(dd2);
  m2(dd2) = m2(dd2)-1;
  val = val + a*b*integrate_data( data, dx_vec, corner_vec, size_vec, m2, n, mask);
  
  n2 = n;       %copy n
  a  = n2(dd1); %save what we have to reduce
  n2(dd1) = n2(dd1)-1;
  b = n2(dd2);
  n2(dd2) = n2(dd2)-1;
  val = val + a*b*integrate_data( data, dx_vec, corner_vec, size_vec, m, n2, mask);
  
  m2 = m;
  n2 = n;
  a  = m2(dd1); %save what we have to reduce
  m2(dd1) = m2(dd1)-1;
  b = n2(dd2);
  n2(dd2) = n2(dd2)-1;
  val = val - a*b*integrate_data( data, dx_vec, corner_vec, size_vec, m2, n2, mask);
  
  m2 = m;
  n2 = n;
  a  = m2(dd2); %save what we have to reduce
  m2(dd2) = m2(dd2)-1;
  b = n2(dd1);
  n2(dd1) = n2(dd1)-1;
  val = val - a*b*integrate_data( data, dx_vec, corner_vec, size_vec, m2, n2, mask);
  
  
  %Now we can consider mixed derivatives. One on mask and one on
  %polynomial.
  m2 = m;      %copy m
  a  = m2(dd1); %save what we have to reduce
  m2(dd1) = m2(dd1)-1;
  val = val + a*integrate_data( data, dx_vec, corner_vec, size_vec, m2, n, mask_gradient{dd2});
  
  n2 = n;      %copy n
  a  = n2(dd1); %save what we have to reduce
  n2(dd1) = n2(dd1)-1;
  val = val - a*integrate_data( data, dx_vec, corner_vec, size_vec, m, n2, mask_gradient{dd2});
  
  m2 = m;      %copy m
  a  = m2(dd2); %save what we have to reduce
  m2(dd2) = m2(dd2)-1;
  val = val + a*integrate_data( data, dx_vec, corner_vec, size_vec, m2, n, mask_gradient{dd1});
  
  n2 = n;      %copy n
  a  = n2(dd2); %save what we have to reduce
  n2(dd2) = n2(dd2)-1;
  val = val - a*integrate_data( data, dx_vec, corner_vec, size_vec, m, n2, mask_gradient{dd1});
end