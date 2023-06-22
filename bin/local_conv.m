function weight = local_conv( func, q, par, gp_vec )

  weight = 0*func;

  dx = par(1);
  dt = par(3);
  
  gpx = gp_vec(1);
  gpy = gp_vec(2);
  gpt = gp_vec(3);

  for i = 1:gpx
	  i %Just had this to track the speed of processing
  for j = 1:gpy
  for k = 1:gpt
    for a = (i + (-q:q))
      if(a<1 || a>gpx)
        continue
      end
    for b = (j + (-q:q))
      if(b<1 || b>gpy)
        continue
      end
    for c = (k + (-q:q))
      if(c<1 || c>gpt)
        continue
      end
      %[a,b,c]
      del_x = (i - a)*dx;
      del_y = (j - b)*dx;
      del_t = (k - c)*dt;

      weight(i,j,k) = weight(i,j,k) + ...
	                exp(-del_x.*del_x/par(1) - del_y.*del_y/par(2) -del_t.*del_t/par(3) )*func(a,b,c);
    end %c
    end %b
    end %a
  end  
  end
  end
end
