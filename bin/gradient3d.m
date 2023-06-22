function [u_y, u_x, u_t] = gradient3d( u, dx_vec )
  %Pass in a scalar with 3 dimensions
  size_vector = size(u);
  
  u_x  = zeros(size_vector);
  u_y  = zeros(size_vector);
  u_t  = zeros(size_vector);

  %https://www.mathworks.com/matlabcentral/answers/16996-assign-multiple-variables  
  [dy dx dt] = feval(@(x) x{:}, num2cell(dx_vec));

  dim = size(u);
  for j = 1:dim(1)
    for i = 1:dim(2)
      for k = 1:dim(3)

        if i == 1
          %u_x(j,i,k) = (-3*u(j,1,k) + 4*u(j,2,k) - u(j,3,k))/(2*dx);
          u_x(j,i,k) = (u(j,i+1,k) - u(j,i,k))/(dx);
        elseif i == dim(2)
          %u_x(j,i,k) = (3*u(j,dim(2),k) - 4*u(j,dim(2)-1,k) + u(j,dim(2)-2,k))/(2*dx);
          u_x(j,i,k) = (-u(j,i-1,k) + u(j,i,k))/(dx);
        else
          u_x(j,i,k) = (u(j,i+1,k) - u(j,i-1,k))/(2*dx);
          %u_x(j,i,k) = (u(j,i+1,k) - u(j,i,k))/(dx);
        end

        if j == 1
          %u_y(j,i,k) = (-3*u(1,i,k) + 4*u(2,i,k) - u(3,i,k))/(2*dy);
          u_y(j,i,k) = (u(j+1,i,k) - u(j,i,k))/(dy);
        elseif j == dim(1)
          %u_y(j,i,k) = (3*u(dim(1),i,k) - 4*u(dim(1)-1,i,k) + u(dim(1)-2,i,k))/(2*dy);
          u_y(j,i,k) = (-u(j-1,i,k) + u(j,i,k))/(dy);
        else
          u_y(j,i,k) = ( u(j+1,i,k) - u(j-1,i,k) )/(2*dy);
          %u_y(j,i,k) = (u(j+1,i,k) - u(j,i,k))/(dy);
        end

	    if k == 1
          u_t(j,i,k) = (u(j,i,k+1) - u(j,i,k))/(dt);
          %u_t(j,i,k) = (-3*u(j,i,1) + 4*u(j,i,2) - u(j,i,3))/(2*dt);
        elseif k == dim(3)
          %u_t(j,i,k) = (3*u(j,i,dim(3)) - 4*u(j,i,dim(3)-1) + u(j,i,dim(3)-2))/(2*dt);
          u_t(j,i,k) = (-u(j,i,k-1) + u(j,i,k))/(dt);
        else
          u_t(j,i,k) = ( u(j,i,k+1) - u(j,i,k-1) )/(2*dt);
          %u_t(j,i,k) = (u(j,i,k+1) - u(j,i,k))/(dt);
        end
      end %k
    end %i
  end %j
end %func