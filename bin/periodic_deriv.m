function d = periodic_deriv( phi, dim, dx)
  d = circshift(phi, -1, dim) - circshift(phi, 1, dim);
  d = d/(2*dx);
end
