function [nx_x, nx_y, nx_t, ny_x, ny_y, ny_t] = nematic_derivative( nx, ny, dx_vec )
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
  nx_x = 0*nx;
  nx_y = 0*nx;
  ny_x = 0*ny;
  ny_y = 0*ny;
  nx_t = 0*ny;
  ny_t = 0*ny;

  %https://www.mathworks.com/matlabcentral/answers/16996-assign-multiple-variables  
  [dy, dx, dt] = feval(@(x) x{:}, num2cell(dx_vec));

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
        past   = [nx(j,i,k-1) ny(j,i,k-1)];
        future = [nx(j,i,k+1) ny(j,i,k+1)];
        
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
        if center*past' < 0
          past = -past;
        end
        if center*future' < 0
          future = -future;
        end
        
        nx_x(j,i,k) = (top(1)   - bottom(1))/(2*dx);
        nx_y(j,i,k) = (right(1) - left(1)  )/(2*dx);
        nx_t(j,i,k) = (future(1) - past(1) )/(2*dt);
        
        ny_x(j,i,k) = (top(2)   - bottom(2))/(2*dx);
        ny_y(j,i,k) = (right(2) - left(2)  )/(2*dx);
        ny_t(j,i,k) = (future(2) - past(2) )/(2*dt);
      end %k
    end %i
  end %j
  
  %Orthogonalize with n
  temp = nx.*nx_x + ny.*ny_x;
  nx_x = nx_x - temp.*nx;
  ny_x = ny_x - temp.*ny;
  
  
  temp = nx.*nx_y + ny.*ny_y;
  nx_y = nx_y - temp.*nx;
  ny_y = ny_y - temp.*ny;
  
  
  temp = nx.*nx_t + ny.*ny_t;
  nx_t = nx_t - temp.*nx;
  ny_t = ny_t - temp.*ny;
end %func
