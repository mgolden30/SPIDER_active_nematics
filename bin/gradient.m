function [u_x, u_y] = gradient( u, dx )
  % IN:
  % u  - a 2D scalar on a uniform grid
  % dx - grid spacing 
  %
  % OUT:
  % u_x - x derivative of u
  % u_y - y derivative of u
  u_x  = 0*u;
  u_y  = 0*u;

  dim = size(u);
  for j = 1:dim(1)
    for i = 1:dim(2)

      if i == 1
        u_x(j,i) = (-3*u(j,1) + 4*u(j,2) - u(j,3))/(2*dx);
      elseif i == dim(2)
        u_x(j,i) = (3*u(j,dim(2)) - 4*u(j,dim(2)-1) + u(j,dim(2)-2))/(2*dx);
      else
        u_x(j,i) = (u(j,i+1) - u(j,i-1))/(2*dx);
      end

      if j == 1
        u_y(j,i) = (-3*u(1,i) + 4*u(2,i) - u(3,i))/(2*dx);
      elseif j == dim(1)
        u_y(j,i) = (3*u(dim(1),i) - 4*u(dim(1)-1,i) + u(dim(1)-2,i))/(2*dx);
      else
        u_y(j,i) = ( u(j+1,i) - u(j-1,i) )/(2*dx);
      end

    end%i
  end%j
end%function
