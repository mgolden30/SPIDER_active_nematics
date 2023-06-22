function [nx_lap, ny_lap] = nematic_laplacian( nx, ny, dx_vec )
  %{
  PURPOSE: Take a consitent derivative of a director field using finite
  differencing. This is done by flipping n in the stencil so that it aligns
  with the center vector.

  INPUT:
  nx - x component of director field
  ny - y component of director field
  dx_vec - [dy dx dt]

  OUTPUT:
  derivatives of n

  %}
  nx_lap = 0*nx;
  ny_lap = 0*ny;
  
  %https://www.mathworks.com/matlabcentral/answers/16996-assign-multiple-variables  
  [dy dx dt] = feval(@(x) x{:}, num2cell(dx_vec));

  dim = size(nx);
  for j = 2:dim(1)-1
    for i = 2:dim(2)-1
      for k = 2:dim(3)-1
        %copy a stencil from nx and ny
        center = [nx(j,i,k) ny(j,i,k)];
        top    = [nx(j,i+1,k) ny(j,i+1,k)];
        bottom = [nx(j,i-1,k) ny(j,i-1,k)];
        left   = [nx(j-1,i,k) ny(j-1,i,k)];
        right  = [nx(j+1,i,k) ny(j+1,i,k)];
        
        %Now we will flip vectors as needed.
        %center*top' is a dot product.
        if center*top' < 0
          top = -top;
        end
        if center*bottom' < 0
          bottom = -bottom;
        end
        if center*right' < 0
          right = -right;
        end
        if center*left' < 0
          left = -left;
        end
          
        nx_lap(j,i,k) = ( top(1) + bottom(1) + left(1) + right(1) - 4*center(1) )/(dx*dx);
        ny_lap(j,i,k) = ( top(2) + bottom(2) + left(2) + right(2) - 4*center(2) )/(dx*dx);
      end %k
    end %i
  end %j
end %func
