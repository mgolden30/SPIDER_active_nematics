function [phi, psi] = nhhd(u, v, dx)
  %Computes the natural 2D Helmholtz-Hodge decomposition
  [gpx, gpy] = size(u);
  
  [div, curl] = div_curl(u,v,dx);
  phi = 0*u;
  psi = 0*u;

  [gpx, gpy] = size(phi);
  for i = 1:gpx
    for j = 1:gpy
      for ii = 1:gpx %integration variable
	for jj = 1:gpy %integration variable
          r = dx*sqrt((i-ii)*(i-ii) + (j-jj)*(j-jj));
	  if(r == 0)
            continue;
	  end
          g_inf = 1/(2*pi) * log(r);
          phi(j,i) = phi(j,i) + dx*dx*g_inf*div(jj,ii);
          psi(j,i) = psi(j,i) - dx*dx*g_inf*curl(jj,ii);
        end
      end
    end
  end

  %Add a constant so that the integral is zero
  average_phi = sum(phi, 'all')/(gpx*gpy);
  average_psi = sum(psi, 'all')/(gpx*gpy);

  phi = phi - average_phi*ones(gpx, gpy);
  psi = psi - average_psi*ones(gpx, gpy);
end
