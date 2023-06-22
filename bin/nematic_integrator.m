function [Q_xx_timeseries, Q_xy_timeseries] = nematic_integrator( u, v, Q_xx_0, Q_xy_0, dx, dt, gpt, out_every ) 
  %This function takes some Q_{ij} and integrates it in time
  %with respect to a background flow
  
  [gpx, gpy] = size(Q_xx_0);

  k1_xx = zeros(gpx,gpy);
  k2_xx = zeros(gpx,gpy);
  k1_xy = zeros(gpx,gpy);
  k2_xy = zeros(gpx,gpy);


  Qrk1_xx = Q_xx_0;
  Qrk2_xx = Q_xx_0;
  Qrk1_xy = Q_xy_0;
  Qrk2_xy = Q_xy_0;
 
  Q_xx_timeseries = zeros(gpx,gpy,gpt);
  Q_xy_timeseries = zeros(gpx,gpy,gpt);
  Q_xx_timeseries(:,:,1) = Q_xx_0;
  Q_xy_timeseries(:,:,1) = Q_xy_0;

  magnitude_n = [];

  for timesteps = 2:gpt;
    %Rk2 is quick and easy
    [k1_xx, k1_xy] = nematic_time_deriv(Qrk1_xx, Qrk1_xy, u, v);
    Qrk2_xx = Qrk1_xx + dt*k1_xx;
    Qrk2_xy = Qrk1_xy + dt*k1_xy;

    [k2_xx, k2_xy] = nematic_time_deriv(Qrk2_xx, Qrk2_xy, u, v);
    Qrk1_xx = Qrk1_xx + dt/2*(k1_xx + k2_xx);
    Qrk1_xy = Qrk1_xy + dt/2*(k1_xy + k2_xy);

    [Qrk1_xx, Qrk1_xy] = nematic_normalize(Qrk1_xx, Qrk1_xy);


    if( mod(timesteps-2,out_every) == 0 )
      plot_nematic( Qrk1_xx, Qrk1_xy, k1_xx, k1_xy, timesteps );
    end

    Q_xx_timeseries(:,:,timesteps) = Qrk1_xx;
    Q_xy_timeseries(:,:,timesteps) = Qrk1_xy;



    %Make sure 
    const = 4*(Qrk1_xx.*Qrk1_xx + Qrk1_xy.*Qrk1_xy);
    next_mag = mean( const, [1 2] );
    magnitude_n = [magnitude_n next_mag];
  end

  plot_constraint( magnitude_n );
end
