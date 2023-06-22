function new_state = dealias(state)
  [points, ~, ~] = size(state);
  u    = state(:,:,1);
  v    = state(:,:,2);
  Q_11 = state(:,:,3);
  Q_12 = state(:,:,4);

  up   = (points/2-1);
  down = -points/2;
  [ky, kx] = meshgrid( down:up, down:up );

  fft_u    = fftshift( fft2(u   ) );
  fft_v    = fftshift( fft2(v   ) );
  fft_Q_11 = fftshift( fft2(Q_11) );
  fft_Q_12 = fftshift( fft2(Q_12) );

  %kill all of the mode above 1/2 the max wavenumber of points/2
  killx = abs(kx) > points/4;
  killy = abs(ky) > points/4;
  kill = (killx + killy ~= 0 );

  fft_u(kill)    = 0;
  fft_v(kill)    = 0;
  fft_Q_11(kill) = 0;
  fft_Q_12(kill) = 0;
  
  
  u    = real( ifft2( ifftshift( fft_u    ) ) );
  v    = real( ifft2( ifftshift( fft_v    ) ) );
  Q_11 = real( ifft2( ifftshift( fft_Q_11 ) ) );
  Q_12 = real( ifft2( ifftshift( fft_Q_12 ) ) );
  
  new_state = state;
  new_state(:,:,1) = u;
  new_state(:,:,2) = v;
  new_state(:,:,3) = Q_11;    
  new_state(:,:,4) = Q_12;
end