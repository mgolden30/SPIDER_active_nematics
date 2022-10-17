addpath('../bin')

clear functions

fprintf('Starting antisymmetric tensor regression...\n');

num_lib     = 9;
num_windows = 10*num_lib;
G = zeros(num_windows, num_lib);

labels = cell(num_lib, 1);
%Loop over windows and integrate.

m = 4*[1 1 1];
n = m;

corners = zeros(3, num_windows);
size_vec= round( 2./dx_vec );
for i=1:num_windows
    edge = 10;
    corners(:,i) = [ randi(Ly-size_vec(1)-2*edge) randi(Lx-size_vec(2)-2*edge) randi(Lt-size_vec(3)-2*edge) ] + [1 1 1]*edge;
end

a = 1;

S  = n1_x + n2_y;
B1 = -(n1.*n1_x + n2.*n1_y);
B2 = -(n1.*n2_x + n2.*n2_y);
un = U.*n1 + V.*n2;

Omega_12 = (U_y - V_x)/2;

S_x = n1_xx + n2_xy; S_y = n1_xy + n2_yy;
[B1_x, B1_y, ~] = gradient3d(B1, dx_vec);
[B2_x, B2_y, ~] = gradient3d(B2, dx_vec);





G_full = [];
lengths = size_vec.*dx_vec;
weights = { 0*x+1, cos( pi*x/lengths(2) ), cos( pi*y/lengths(1)), cos( pi*t/lengths(3)) };
mask0 = mask;

for ii = 1:numel(weights)
a = 1;

clear functions            %We will recompute derivatives of the mask
mask = mask0.*weights{ii}; %Multiply your weight by anything you want



s_12 = n1.*n2_t; 
s_21 = n2.*n1_t;
G(:,a) = gpu_integrate_antisymm_tensor( s_12, s_21, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_{[i} \partial_t n_{j]}";
a=a+1;

s_12 = S.*U.*n2; 
s_21 = S.*V.*n1;
G(:,a) = gpu_integrate_antisymm_tensor( s_12, s_21, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S u_{[i} n_{j]}";
a=a+1;

s_12 = n1.*S_y; 
s_21 = n2.*S_x;
G(:,a) = gpu_integrate_antisymm_tensor( s_12, s_21, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_{[i} \nabla_{j]} S";
a=a+1;

s_12 = S.*n1.*B2; 
s_21 = S.*n2.*B1;
G(:,a) = gpu_integrate_antisymm_tensor( s_12, s_21, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S n_{[i} B_{j]}";
a=a+1;

s_12 = n1.*n1.*U_y + n2.*n1.*V_y; 
s_21 = n1.*n2.*U_x + n2.*n2.*V_x;
G(:,a) = gpu_integrate_antisymm_tensor( s_12, s_21, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_k n_{[i} \nabla_{j]} u_k";
a=a+1;

s_12 = Omega_12; 
s_21 = -Omega_12;
G(:,a) = gpu_integrate_antisymm_tensor( s_12, s_21, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\Omega_{ij}";
a=a+1;

%{
s_12 = un.*B1.*n2; 
s_21 = un.*B2.*n1;
G(:,a) = gpu_integrate_antisymm_tensor( s_12, s_21, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_k n_k B_{[i} n_{j]}";
a=a+1;
%}

s_12 = U.*B2; 
s_21 = V.*B1;
G(:,a) = gpu_integrate_antisymm_tensor( s_12, s_21, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_{[i} B_{j]}";
a=a+1;

s_12 = n1.*( n1.*B2_x + n2.*B2_y );
s_21 = n2.*( n1.*B1_x + n2.*B1_y );
G(:,a) = gpu_integrate_antisymm_tensor( s_12, s_21, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_{[i|} n_k \nabla_k B_{|j]}";
a=a+1;


s_12 = B2_x; 
s_21 = B1_y;
G(:,a) = gpu_integrate_antisymm_tensor( s_12, s_21, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\nabla_{[i} B_{j]}";
a=a+1;


G_full = [G_full; G];
end %end weight iteration

G = G_full;
mask = mask0;
fprintf("done integrating.\n")









function [s_11, s_12, s_21, s_22] = exterior_product(n1,n2,m1,m2)
  %symmetrized exterior product
  s_11 = n1.*m1;
  s_12 = (n1.*m2 + n2.*m1)/2;
  s_21 = s_12;
  s_22 = n2.*m2;
end


function vals = gpu_integrate_antisymm_tensor( s_12, s_21, dx_vec, corners, size_vec, m, n, mask )
  %Only two degrees of freedom
  s_12 = (s_12 - s_21)/2; %compute antisymmetric part
  vals = gpu_integrate_data( s_12, [], dx_vec, corners, size_vec, m, n, mask );
  vals = gather(vals);
end