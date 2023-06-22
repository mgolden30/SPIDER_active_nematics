function smoother_data = moving_least_squares_2d( data, n)
  %{
  PURPOSE:
  I hate noisy data. This function performs a local quadratic fit to the
  data.

  INPUT:
  data - 2D data
  n - number of neighbors on each side to use. For example, n=3 will use 3
  neighbors on each side resulting in a 7x7 box of samples for the
  polynomial fit. Something this order is ideal as you want to avoid
  overfitting.
  

  OUTPUT:
  smoother_data - exactly what it says it is
  %}

  [a,b] = size(data);
  smoother_data = zeros( size(data) );
  
  parfor i=1:a
    for j=1:b
      stencil_x = i + (-n:n); stencil_x( stencil_x < 1 | stencil_x > a) = [];
      stencil_y = j + (-n:n); stencil_y( stencil_y < 1 | stencil_y > b) = [];
      
      [sx, sy] = meshgrid( stencil_x - i, stencil_y - j ); %shift point of interest to the origin
      sx = reshape( sx, [numel(sx),1] );
      sy = reshape( sy, [numel(sy),1] );

      local_data = data(stencil_x, stencil_y);
      local_data = reshape( local_data, [numel(local_data),1] );
      
      A = [ 0*sx+1, sx, sy, sx.^2, sy.^2, sx.*sy ]; %All of our basis functions go here
      coeff = linsolve( A'*A, A'*local_data);
      smoother_data(i,j) = coeff(1); %since our point of interest is the origin, just read out the constant coefficient.
    end
  end
end