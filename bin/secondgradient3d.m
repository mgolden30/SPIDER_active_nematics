function [u_yy, u_xx, u_tt, u_xy, u_xt, u_yt] = secondgradient3d( u, dx_vec )
  %Pass in a scalar with 3 dimensions
  size_vector = size(u);
  
  u_xx = zeros(size_vector);
  u_yy = zeros(size_vector);
  u_tt = zeros(size_vector);
  
  u_xt = zeros(size_vector);
  u_yt = zeros(size_vector);
  u_xy = zeros(size_vector);
  
  %https://www.mathworks.com/matlabcentral/answers/16996-assign-multiple-variables  
  [dy dx dt] = feval(@(x) x{:}, num2cell(dx_vec));

  dim = size(u);
  for j = 2:dim(2)-1
    for i = 2:dim(1)-1
      for k = 2:dim(3)-1
        u_xx(j,i,k) = ( u(j,i+1,k) - 2*u(j,i,k) + u(j,i-1,k) )/(dx*dx);
        u_yy(j,i,k) = ( u(j+1,i,k) - 2*u(j,i,k) + u(j-1,i,k) )/(dx*dx);
        u_tt(j,i,k) = ( u(j,i,k+1) - 2*u(j,i,k) + u(j,i,k+1) )/(dt*dt);
        
        u_xy(j,i,k) = ( u(j+1,i+1,k) + u(j-1,i-1,k) - u(j-1,i+1,k) - u(j+1,i-1,k) )/(4*dx*dx);
        u_xt(j,i,k) = ( u(j,i+1,k+1) + u(j,i-1,k-1) - u(j,i+1,k-1) - u(j,i-1,k+1) )/(4*dx*dt);
        u_yt(j,i,k) = ( u(j+1,i,k+1) + u(j-1,i,k-1) - u(j+1,i,k-1) - u(j-1,i,k+1) )/(4*dx*dt);
      end %k
    end %i
  end %j
end %func