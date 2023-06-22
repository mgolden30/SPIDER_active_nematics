function val = integrate_derivs( data, deriv_vec, dx_vec, corner_vec, size_vec, weight_obj )
  %{
   PURPOSE:
   This function was created to save me some development time. I am using
   my weight class I wrote to generalize integration by parts without
   hurting myself.
  
   This function evaluates the integral:
       weight*(bunch of derivatives)data
  %}

  [dy, dx, dt] = feval(@(x) x{:}, num2cell(dx_vec));
  [~,N] = size(deriv_vec);
  %recursively apply the derivatives
  temp = weight_obj;
  for i = 1:N
    temp = temp.derivative( deriv_vec(i) );
  end
  
  
  %Some shorthand for the components of corner_vec
  [cy, cx, ct] = feval(@(x) x{:}, num2cell(corner_vec));
  
  %restrict data to the region of integration
  data = data( cy:(cy+size_vec(1)), cx:(cx+size_vec(2)), ct:(ct+size_vec(3)));
  
  %evaluate the weight object in this domain
  w    = temp.eval(corner_vec, size_vec);
  
  val = sum( data.*w, 'all' );
  val = val*((-1)^N)*dx*dy*dt; %Lastly correct for the change in sign from integration by parts
end