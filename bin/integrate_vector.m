function val = integrate_vector( u, v, dx, corner_vec, size_vec, m_column, n_column, mask )
  %{
  so we need a weight vector <w_y, -w_x> that is the curl of a scalar w. Then
  \int dA <u,v> \cdot <w_y, -w_x>  = \int dA (v_x - u_y) w  if w and its derivatives vanish on the boundary.
  %}
  val =       integrate_deriv( v, 2, dx, corner_vec, size_vec, m_column, n_column, mask );
  val = val - integrate_deriv( u, 1, dx, corner_vec, size_vec, m_column, n_column, mask );
end
