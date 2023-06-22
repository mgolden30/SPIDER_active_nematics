function smoothed_data = smooth_out( data )
  [ny, nx, nt] = size(data);
  smoothed_data = data;

  for i = 2:(ny-1)
    for j = 2:(nx-1)
      for k = 2:(nt-1)
        smoothed_data(i,j,k) = (data(i,j,k) + ...
	                       data(i+1,j,k) + data(i-1,j,k) + ...
	                       data(i,j+1,k) + data(i,j-1,k) + ...
	                       data(i,j,k+1) + data(i,j,k-1))/7;
      end
    end
  end


end
