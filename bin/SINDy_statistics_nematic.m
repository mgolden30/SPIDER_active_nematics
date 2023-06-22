function [xi_average, xi_stdev, eta_average, eta_stdev, gamma_vec ] = SINDy_statistics_nematic( theta, dXdt )
  num_windows     = 100;
  num_lib         = 8; %Usually 8
  nw = num_windows;
  nl = num_lib;
  xi_average  = [];
  xi_stdev    = [];
  eta_average = [];
  eta_stdev   = [];
  gamma_vec   = [];

  %Note 10 ^ -0.3 ~ 0.5
  for log_gamma = -3:0.1:-0.3 %base 10
    gamma = 10^log_gamma;
    subsamples = 100;
    xis  = cell(1, subsamples);
    etas = cell(1, subsamples);

    for i=1:subsamples
      y = randsample(nw, 50);
      small_dXdt  = zeros(1, 2*50);
      small_theta = zeros(2*50, num_lib);

      small_dXdt(y)      =dXdt(y);
      small_dXdt(y+50)   =dXdt(y+nw);
      small_theta(y,:)   =theta(y,:);
      small_theta(y+50,:)=theta(y+nw,:);

      size(small_theta)
      size(small_dXdt)
      [xis{i}, etas{i}]= SINDy_var_gamma( small_theta, small_dXdt', gamma );
    end

    %compute the average
    xi_av = zeros(num_lib);
    eta_av = 0;
    for i=1:subsamples
      xi_av  = xi_av + xis{i};
      eta_av = eta_av + etas{i};
    end
    xi_av  = xi_av/subsamples;
    eta_av = eta_av/subsamples;

    %compute std_dev
    xi_sd  = zeros(num_lib);
    eta_sd = 0;
    for i=1:subsamples
      xi_sd  =  xi_sd + (xis{i}  -  xi_av).^2;
      eta_sd = eta_sd + (etas{i} - eta_av).^2;
    end
    xi_sd  = sqrt( xi_sd/(subsamples-1));
    eta_sd = sqrt(eta_sd/(subsamples-1));

    xi_average = [xi_average xi_av];
    xi_stdev   = [xi_stdev xi_sd];
    eta_average = [eta_average eta_av];
    eta_stdev   = [eta_stdev   eta_sd];
    gamma_vec   = [gamma_vec   gamma];
  end %log_gamma loop

end

