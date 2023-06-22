function val = integrate_data( data, dx_vec, corner_vec, size_vec, m, n, mask )
  % PURPOSE:
  % The purpose of this function is to integrate scalar data over a spacetime cube with some weight function
  %
  % IN:
  % data   - some spatio-temporal scalar data.
  % dx_vec - a vector of displacements in each direction
  % corner_vec - the bottom corner of some spacetime cube (in gridpoints)
  % size_vec - the dimensions of the spacetime cube (in gridpoints)
  % m - powers of the polynomial weight
  % n - other powers of the polynomial weight
  % mask - another weight function that can be manually input
  %
  % OUT:
  % val - the integrated value
  
  [dy, dx, dt] = feval(@(x) x{:}, num2cell(dx_vec));
  [cy, cx, ct] = feval(@(x) x{:}, num2cell(corner_vec));
  [my, mx, mt] = feval(@(x) x{:}, num2cell(m));
  [ny, nx, nt] = feval(@(x) x{:}, num2cell(n));
  [sy, sx, sz] = size(data);
  
  %first restict stuff to the region we care about
  data = data( cy:(cy+size_vec(1)), cx:(cx+size_vec(2)), ct:(ct+size_vec(3)));
  mask = mask( cy:(cy+size_vec(1)), cx:(cx+size_vec(2)), ct:(ct+size_vec(3)));
  
  to_be_integrated = data.*mask;
  [X, Y, T] = meshgrid( 0:size_vec(1), 0:size_vec(2), 0:size_vec(3));

  Y = Y*dy;
  X = X*dx;
  T = T*dt;
  
  polynomial_weight = (Y.^my).*(size_vec(1)*dy - Y).^ny.*...
                      (X.^mx).*(size_vec(2)*dx - X).^nx.*...
                      (T.^mt).*(size_vec(3)*dt - T).^nt;

  to_be_integrated = to_be_integrated.*polynomial_weight;

  %Just sum over all the points now
  val = sum(to_be_integrated, 'all')*dx*dy*dt;
end
