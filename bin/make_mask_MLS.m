function mask = make_mask_MLS( theta, dx )
  %{
  PURPOSE:
  Creates a function "mask" which is 1 far from defects and 0 (exactly)
  near them.
 
  IN:
  Q_11 - component of the nematic tensor
  Q_12 - "
  dx - grid spacing
  alpha and beta - parameters to finely tune the mask.
                   I designed them so that alpha=beta=1 should be a good
                   default choice.

  boolean_mask_branch_cuts - a boolean (true/false) that determines if the
  branch cuts where n changes sign should be masked.

  OUT:
  mask - a 2+1D matrix
  %}
  
  n1 = cos(theta);
  n2 = sin(theta);

  %Calculate gradient
  dx_vec = [dx dx dx];
  
  [n1_x, n1_y, ~, n2_x, n2_y, ~] = nematic_derivative(n1, n2, dx_vec);
  
  %Curvature
  curv = n1_x.^2 + n1_y.^2 + n2_x.^2 + n2_y.^2;
  
  mask = curv < (0.1)^2/dx/dx;
  
  mask = 1-mask;
  SE = strel('disk', 2, 4); %make a structural element
  for i=1:size(mask,3)
    mask(:,:,i) =  imerode(mask(:,:,i), SE);
    for j=1:5
      mask(:,:,i) =  imdilate(mask(:,:,i), SE);
    end
  end
  mask = 1-mask;
  
  
  for t=1:100
    imagesc(mask(:,:,t));
    title(""+t)
    pause(1e-1);
  end

  for i=1:3
    mask( mask < 0.1 ) = 0; %I want exactly 0.
    tic
    mask = moving_least_squares_3d( mask, 2 );
    toc
  end
end