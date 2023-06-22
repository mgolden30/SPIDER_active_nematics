function [u_yy, u_xx, u_tt] = double_gradient3d( u, dx_vec )
  %Pass in a scalar with 3 dimensions
  u_xx  = 0*u;
  u_yy  = 0*u;
  u_tt  = 0*u;

  %https://www.mathworks.com/matlabcentral/answers/16996-assign-multiple-variables  
  [dy dx dt] = feval(@(x) x{:}, num2cell(dx_vec));

  %https://www.dam.brown.edu/people/alcyew/handouts/numdiff.pdf
  dim = size(u);
  for j = 1:dim(1)
    for i = 1:dim(2)
      for k = 1:dim(3)

        if i == 1
          u_xx(j,i,k) = ( 2*u(j,1,k) - 5*u(j,2,k) + 4*u(j,3,k) - u(j,4,k) )/(dx*dx);
        elseif i == dim(2)
          u_xx(j,i,k) = ( 2*u(j,dim(2),k) - 5*u(j,dim(2)-1,k) + 4*u(j,dim(2)-2,k) - u(j,dim(2)-3,k) )/(dx*dx);
        else
          u_xx(j,i,k) = (u(j,i+1,k) - 2*u(j,i,k)  + u(j,i-1,k))/(dx*dx);
        end

        if j == 1
          u_yy(j,i,k) = ( 2*u(1,i,k) - 5*u(2,i,k) + 4*u(3,i,k) - u(4,i,k) )/(dy*dy);
        elseif j == dim(1)
          u_yy(j,i,k) = ( 2*u(dim(1),i,k) - 5*u(dim(1)-1,i,k) + 4*u(dim(1)-2,i,k) - u(dim(1)-3,i,k) )/(dy*dy);
        else
          u_yy(j,i,k) = ( u(j+1,i,k) - 2*u(j,i,k) + u(j-1,i,k) )/(dy*dy);
        end

	    if k == 1
          u_tt(j,i,k) = ( 2*u(j,i,1) - 5*u(j,i,2) + 4*u(j,i,3) - u(j,i,4) )/(dt*dt);
        elseif k == dim(3)
          u_tt(j,i,k) = ( 2*u(j,i,dim(3)) - 5*u(j,i,dim(3)-1) + 4*u(j,i,dim(3)-2) - u(j,i,dim(3)-3) )/(dt*dt);
        else
          u_tt(j,i,k) = ( u(j,i,k+1) - 2*u(i,j,k) + u(j,i,k-1) )/(dt*dt);
        end

      end %k
    end %i
  end %j
end %func