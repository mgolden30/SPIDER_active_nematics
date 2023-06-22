function dfdx = nifty_deriv( f )
  [~, points] = size(f);
  
  dx = 2*pi/points;
  x = 1:points; %eventually replace this with meshgrid.
  x = dx*(x-1);

  g = exp(-1./x - 1./(2*pi - x));

  %This is needed for even points
  up   = (points/2-1);
  down = -points/2;
  kx = down:up;
  kx(1) = 0.; %Set the zig-zag mode to zero

  fg_fft = fftshift(fft(f.*g));
  g_fft  = fftshift(fft(g));
  
  fg_deriv =  fg_fft.*kx*complex(0,1);
   g_deriv =   g_fft.*kx*complex(0,1);
  
  fg_deriv = ifft(ifftshift(fg_deriv));
   g_deriv = ifft(ifftshift( g_deriv));
  
  dfdx = ( fg_deriv - f.*g_deriv )./g;
end
