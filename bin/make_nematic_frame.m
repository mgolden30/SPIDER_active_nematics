function make_nematic_frame( Q_11, Q_12, ii  )
  theta = atan2(Q_12, Q_11)/2;

  nx = cos(theta);
  ny = sin(theta);

  [gpy, gpx] = size(Q_11)

  [x_start, y_start] = meshgrid(1:5:gpx, 1:5:gpy);

  %Turn these into vectors
  [a,b] = size(x_start);
  num_lines = a*b;

  x_start = reshape( x_start, [num_lines 1] );
  y_start = reshape( y_start, [num_lines 1] );

  fibers_x = cell(1, num_lines);
  fibers_y = cell(1, num_lines);

  ffibers_x = cell(1, num_lines);
  ffibers_y = cell(1, num_lines);



  for i=1:num_lines
    x = x_start(i);
    y = y_start(i);

    x_series = [];
    y_series = [];

    vx = 1;
    vy = 1;

    dlambda = -0.1;
    for lambda = 1:100
      previous_vx = vx;
      previous_vy = vy;

      rounded_x = round(x);
      rounded_y = round(y);
 
      vx = nx( rounded_x, rounded_y );
      vy = ny( rounded_x, rounded_y );

      dot = vx*previous_vx + vy*previous_vy;
      if ( dot<0 )
          %Flip vx and vy
	  vx = -vx;
	  vy = -vy;
      end

      %step forward in lambda
      x = x + vx*dlambda;
      y = y + vy*dlambda;

      x_series = [x_series x];
      y_series = [y_series y];

      if(round(x) > gpx || round(y) > gpy || round(x) < 1 || round(y) < 1 )
        break
      end
    end %for loop over lambda
    ffibers_x{i} = x_series;
    ffibers_y{i} = y_series;
  end






  for i=1:num_lines
    x = x_start(i);
    y = y_start(i);

    x_series = [];
    y_series = [];

    vx = 1;
    vy = 1;

    dlambda = 0.1;
    for lambda = 1:100
      previous_vx = vx;
      previous_vy = vy;

      rounded_x = round(x);
      rounded_y = round(y);
 
      vx = nx( rounded_x, rounded_y );
      vy = ny( rounded_x, rounded_y );

      dot = vx*previous_vx + vy*previous_vy;
      if ( dot<0 )
          %Flip vx and vy
	  vx = -vx;
	  vy = -vy;
      end

      %step forward in lambda
      x = x + vx*dlambda;
      y = y + vy*dlambda;

      x_series = [x_series x];
      y_series = [y_series y];

      if(round(x) > gpx || round(y) > gpy || round(x) < 1 || round(y) < 1 )
        break
      end
    end %for loop over lambda
    fibers_x{i} = x_series;
    fibers_y{i} = y_series;
  end

  %Now plot all of these lines
  thick = 0.7;
  plot(fibers_x{1},   fibers_y{1}, 'color', 'blue', 'LineWidth', thick );
  plot(ffibers_x{1}, ffibers_y{1}, 'color', 'blue', 'LineWidth', thick );
  hold on
  for i=2:num_lines
    plot(fibers_x{i},   fibers_y{i}, 'color', 'blue', 'LineWidth', thick );
    plot(ffibers_x{i}, ffibers_y{i}, 'color', 'blue', 'LineWidth', thick );
  end
  hold off

  xlim([1 gpx]);
  ylim([1 gpy]);

  %axis(gca, 'tight')
  set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])


  set(gcf,'PaperPositionMode','auto')
  set(gcf,'position',[0,0,100,100])
  set(gca,'position',[0,0,1,1])

  filename = sprintf('reconstruction/fibers_%03d.png', ii);
  saveas(gcf, filename);

end
