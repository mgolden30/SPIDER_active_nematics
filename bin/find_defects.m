function linear_idx = find_defects( Q_11, Q_12, dx )
  %{
  PURPOSE:
  Finds the topological defects of the nematic field. This is done by looking at the determinant of the gradient of
  the director field. This was suggested by Roman
 
  IN:
  Q_11 - component of the nematic tensor
  Q_12 - "
  dx - grid spacing

  OUT:

  %}

  [gpx,gpy,gpt] = size(Q_11);
  [x,y] = meshgrid( 1:gpx, 1:gpy );
  x = x*dx;
  y = y*dx;

  %Find theta
  theta = atan2(Q_12, Q_11)/2;
  nx = cos(theta);
  ny = sin(theta);

  %Calculate gradient
  nx_x = nematic_derivative(nx, dx, 2);
  ny_x = nematic_derivative(ny, dx, 2);
  nx_y = nematic_derivative(nx, dx, 1);
  ny_y = nematic_derivative(ny, dx, 1);
  det1 = nx_x.*ny_y - nx_y.*ny_x;
  det2 = nx_x.*ny_y + nx_y.*ny_x;

  index1 = abs(det1)<abs(det2);
  index2 = abs(det2)<abs(det1);
  determin = index1.*det1 + index2.*det2;

%  defects = abs(determin) > 10^10;  %Find indiced where defects are
%  linear_idx = find(defects);

%  conv = matfile('correct_diff_convolution.mat');
  conv = matfile('roman_weight_2.mat');

%  weight = local_conv( abs(determin), 3, [dx dx dt], [gpx gpy gpt] );
%  weight = local_conv( weight, 3, [dx dx dt], [gpx gpy gpt] );

%  conv = matfile('tanh_convolution.mat');
  weight = conv.weight;  
  
  Q_11_x = derivative( Q_11, dx, 2);
  Q_11_y = derivative( Q_11, dx, 1);
  Q_12_x = derivative( Q_12, dx, 2);
  Q_12_y = derivative( Q_12, dx, 1);
 

 %{
  q_1 = Q_12_x - Q_11_y;
  q_2 = -Q_11_x - Q_12_y;

  qn_dot = nx.*q_1 + ny.*q_2;
  
  q_1 = q_1 - qn_dot.*nx;
  q_2 = q_2 - qn_dot.*ny;


  determin = q_1.^2 + q_2.^2;
  weight = log(abs( determin(:,:,ii) ));
  %}

  %weight = tanh( exp(25) ./ weight ).^2;
  %for i = 1:10
  %  weight = smooth_out( weight );
  %end
  %save('roman_weight_2.mat', 'weight')
  %return;

  numcontours = 30;
  contourf(x,y, advec(:,:,ii).*weight(:,:,ii), numcontours,'edgecolor','none');
  colormap 'jet';
  %caxis([-1 1]);
  colorbar();

  hold on
  quiver( x, y, nx(:,:,ii), ny(:,:,ii) )
  hold off

  title('smoothed tanh(e^{25}/|det|)^2 * advection')
  output = sprintf('~/Workspace/machine_learning/frames/mask_advection_%03d.png', ii);
  saveas( gcf, output );
end
