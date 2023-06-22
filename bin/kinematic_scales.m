function [length_scale, time_scale] = kinematic_scales( U, V, dx_vec)
  %{
  PURPOSE:
  The purpose of this function is to find the length and time scales such
  that the average velocity and vorticity are O(1). I call such units
  kinematic units.
  %}

  %First filter data
  L = 5;
  U = imgaussfilt3(U, L);  
  V = imgaussfilt3(V, L);
  
  %Throw out border data
  L = 2*L;
  U = U( L:end-L, L:end-L, L:end-L );
  V = V( L:end-L, L:end-L, L:end-L );
  
  [U_y, ~,   ~] = gradient3d(U, dx_vec);
  [~,   V_x, ~] = gradient3d(V, dx_vec);

  omega = V_x - U_y;
  vel_mag = sqrt(U.^2 + V.^2);
  
  U = mean(vel_mag, 'all');
  W = mean(abs(omega), 'all');
  
  length_scale = U/W;
  time_scale   = 1/W;
end