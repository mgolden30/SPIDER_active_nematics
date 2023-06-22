function [smoother_data] = moving_least_squares_3d( data, n)
  %{
  PURPOSE:
  I hate noisy data. This function performs a local quadratic fit to noisy
  data in an attempt to smooth it.

  INPUT:
  data - 3D data
  n - number of neighbors on each side to use. For example, n=3 will use 3
  neighbors on each side resulting in a 7x7x7 box of samples for the
  polynomial fit. Something this order is ideal as you want to avoid
  overfitting.
  

  OUTPUT:
  smoother_data - exactly what it says it is
  f_i - derivatives of the data
  %}

  tic
  
  [a,b,c] = size(data);
  
  smoother_data = zeros( size(data) );
  %f_x = zeros( size(data) );
  %f_y = zeros( size(data) );
  %f_z = zeros( size(data) );
  
  parfor i=1:a
        stencil_x = i + (-n:n); stencil_x( stencil_x < 1 | stencil_x > a) = [];
    for j=1:b
        stencil_y = j + (-n:n); stencil_y( stencil_y < 1 | stencil_y > b) = [];
      for k=1:c
        stencil_z = k + (-n:n); stencil_z( stencil_z < 1 | stencil_z > c) = [];
      
        
        [sx, sy, sz] = meshgrid( stencil_x - i, stencil_y - j, stencil_z - k ); %shift point of interest to the origin
        sx = reshape( sx, [numel(sx),1] );
        sy = reshape( sy, [numel(sy),1] );
        sz = reshape( sz, [numel(sy),1] );

        local_data = data(stencil_x, stencil_y, stencil_z);
        local_data = reshape( local_data, [numel(local_data),1] );
      
        %Make a matrix of our polynomial basis functions
        A = [ 0*sx+1, sx, sy, sz, sx.^2, sy.^2, sz.^2, sx.*sy, sx.*sz, sy.*sz ];

        %coeff = A\local_data;
        coeff = linsolve( A'*A, A'*local_data );
        
        smoother_data(i,j,k) = coeff(1); %since our point of interest is the origin, just read out the constant coefficient.
        
        %We get derivatives for free
        %f_x(i,j,k) = coeff(2)/dx_vec(1);
        %f_y(i,j,k) = coeff(3)/dx_vec(2);
        %f_z(i,j,k) = coeff(4)/dx_vec(3);
      end
    end
  end
  
  toc
end