function vals = gpu_integrate_data( data, derivs, dx_vec, corners, size_vec, m, n, mask )
  %{
  PURPOSE:
  I am rewriting some old code in an attempt to speed up integration times.
  I want to accomplish two things with this integration routine:
  1. Cut down on redundancy. I want to call this function a single time per
     library term. Maybe a couple of times if different components of the term
     can be integrated by parts in various ways.
  2. Leverage GPU acceleration. Sums scale logarithmically in parallel,
     instead of linearly in the sequential case.
  
  INPUT:
  data    - a 2+1D matrix containing the data to integrate with the trapezoid
            rule.

  derivs  - [1,1] for instance uses integration by parts to integrate two
            derivatives along the first dimension (could be x or y
            depending on your convention)

  dx_vec  - [dy dx dt] in my convention

  corners - a 3xnum_windows matrix containing corners of your integration
            domains.

  size_vec- [24 24 48] as an example. The size of the cube you want to integrate over.
  
  m,n     - powers of the polynomial weights. Use [4 4 4] for both when in doubt. 
  
  mask    - an arbitrary function you want to additionally multiply the
            weight by. Use this to get rid of questionable data.

  OUTPUT:
  vals - a vector of integrated values.
  %}

  persistent mask_deriv;
  persistent mask_second_deriv;  

  %This way you can pass a column or row vector.
  d = numel(derivs);
  
  
  %No integration by parts is the base case
  if d == 0 %derivs = [];  
    [three, nw] = size(corners);
    assert(three == 3);
  
    %First move data to the gpu.
    gpu_data = gpuArray(data);
    gpu_mask = gpuArray(mask);
  
    vals = zeros(nw,1);
    vals = gpuArray(vals);
         
    %Multiply our data by the mask
    gpu_data = gpu_data.*gpu_mask;
  
    %All coordinates should be from 0 to 1
    idx_x   = 0:size_vec(1);
    idx_y   = 0:size_vec(2);
    idx_t   = 0:size_vec(3);
  
    [x,y,t] = meshgrid( idx_x, idx_y, idx_t );
  
    %move all these guys over to the GPU so that calculations are faster
    x = gpuArray(x);
    y = gpuArray(y);
    t = gpuArray(t);
  
    x = dx_vec(2) * x;
    y = dx_vec(1) * y;
    t = dx_vec(3) * t;
    polynomial_weight = (x.^m(2)).*(size_vec(2)*dx_vec(2) - x).^n(2) .*... 
                        (y.^m(1)).*(size_vec(1)*dx_vec(1) - y).^n(1) .*...
                        (t.^m(3)).*(size_vec(3)*dx_vec(3) - t).^n(3);
  
    for i = 1:nw
      small_data = gpu_data( corners(1,i) + idx_x, corners(2,i) + idx_y, corners(3,i) + idx_t );
      vals(i) = sum( small_data.*polynomial_weight, 'all' );
    end
  
    scaling = prod( dx_vec );
    vals = vals*scaling;
    
    vals = gather(vals);
    return;
  end
  
  if d==1
    if isempty(mask_deriv)
      mask_deriv = cell(3,1);
      fprintf('generating derivatives of mask\n')
      [ mask_deriv{1}, mask_deriv{2}, mask_deriv{3} ] = gradient3d( mask, dx_vec );
    end
    
    m2 = m;
    m2(derivs(1)) = m2(derivs(1)) - 1;
    n2 = n;
    n2(derivs(1)) = n2(derivs(1)) - 1;
    
    vals1 =              -gpu_integrate_data( data, [], dx_vec, corners, size_vec, m,  n, mask_deriv{derivs(1)} );
    vals2 = -m(derivs(1))*gpu_integrate_data( data, [], dx_vec, corners, size_vec, m2, n, mask );
    vals3 =  n(derivs(1))*gpu_integrate_data( data, [], dx_vec, corners, size_vec, m, n2, mask );
    vals = vals1 + vals2 + vals3;
  end
  
  if d==2
    %{
    For two derivatives, there are 9 terms in the product rule.
    %}
    d1 = derivs(1);
    d2 = derivs(2);
      
    if isempty( mask_second_deriv )
      mask_second_deriv = cell(3,3);
      fprintf('generating second derivatives of mask\n')
      [mask_second_deriv{1,1}, mask_second_deriv{2,2}, mask_second_deriv{3,3}, mask_second_deriv{1,2}, mask_second_deriv{2,3}, mask_second_deriv{1,3}]... 
          = secondgradient3d( mask, dx_vec );
      mask_second_deriv{2,1} = mask_second_deriv {1,2};
      mask_second_deriv{3,1} = mask_second_deriv {1,3};
      mask_second_deriv{3,2} = mask_second_deriv {2,3};
    end
    
    if isempty(mask_deriv)
      mask_deriv = cell(3,1);
      fprintf('generating derivatives of mask')
      [ mask_deriv{1}, mask_deriv{2}, mask_deriv{3} ] = gradient3d( mask, dx_vec );
    end
    
    m_1 = m;
    m_1(d1) = m_1(d1) - 1;
    
    n_1 = n;
    n_1(d1) = n_1(d1) - 1;
    
    m_2 = m;
    m_2(d2) = m_2(d2) - 1;
    
    n_2 = n;
    n_2(d2) = n_2(d2) - 1;
    
    m_12 = m_1;
    m_12(d2) = m_12(d2) - 1;
    
    n_12 = n_1;
    n_12(d2) = n_12(d2) - 1;
    
    
    val1 =        gpu_integrate_data( data, [], dx_vec, corners, size_vec, m,    n, mask_second_deriv{d1,d2} );
    
    val2 =  m(d2)*gpu_integrate_data( data, [], dx_vec, corners, size_vec, m_2,  n, mask_deriv{d1} );
    val3 = -n(d2)*gpu_integrate_data( data, [], dx_vec, corners, size_vec, m,  n_2, mask_deriv{d1} );
    
    val4 =  m(d1)*gpu_integrate_data( data, [], dx_vec, corners, size_vec, m_1,  n, mask_deriv{d2} );
    val5 = -n(d1)*gpu_integrate_data( data, [], dx_vec, corners, size_vec, m,  n_1, mask_deriv{d2} );
    
    val6 = -m(d1)*n(d2)*gpu_integrate_data( data, [], dx_vec, corners, size_vec, m_1,  n_2, mask );
    val7 = -m(d2)*n(d1)*gpu_integrate_data( data, [], dx_vec, corners, size_vec, m_2,  n_1, mask );
    
    if d1 == d2
      scale8 = m(d1)*(m(d1)-1);
      scale9 = n(d1)*(n(d1)-1);
    else
      scale8 = m(d1)*m(d2);
      scale9 = n(d1)*n(d2);
    end
    
    val8 = scale8*gpu_integrate_data( data, [], dx_vec, corners, size_vec, m_12,  n, mask );
    val9 = scale9*gpu_integrate_data( data, [], dx_vec, corners, size_vec, m,  n_12, mask );
    
    %val_matrix = [val1 val2 val3 val4 val5 val6 val7 val8 val9 ];
    %vecnorm(val_matrix)
    
    vals = val1 + val2 + val3 + val4 + val5 + val6 + val7 + val8 + val9;
  end
  
  vals = gather(vals);
  
  if d>2
    error('too many derivatives. Not supported yet.') 
  end
end