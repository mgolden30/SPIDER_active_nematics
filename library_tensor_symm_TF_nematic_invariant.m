addpath('../bin')

clear functions

fprintf('Starting symmetric trace-free tensor regression...\n');

num_lib     = 14;
num_windows = 10*num_lib;
G = zeros(2*num_windows, num_lib);

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



G_full = [];
lengths = size_vec.*dx_vec;
weights = { 0*x+1, cos( pi*x/lengths(2) ), cos( pi*y/lengths(1)), cos( pi*t/lengths(3)) };
mask0 = mask;

for ii = 1:numel(weights)

a = 1;

clear functions            %We will recompute derivatives of the mask
mask = mask0.*weights{ii}; %Multiply your weight by anything you want

S  = n1_x + n2_y;
B1 = -(n1.*n1_x + n2.*n1_y);
B2 = -(n1.*n2_x + n2.*n2_y);
un = U.*n1 + V.*n2;

Omega_12 = (U_y - V_x)/2;

S_x = n1_xx + n2_xy; s_y = n1_xy + n2_yy;
[B1_x, B1_y, ~] = gradient3d(B1, dx_vec);
[B2_x, B2_y, ~] = gradient3d(B2, dx_vec);


[s_11, s_12, s_21, s_22] = exterior_product(n1,n2,n1,n2);
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "Q'_{ij}";
a=a+1;

[s_11, s_12, s_21, s_22] = exterior_product(n1,n2,S.*S.*n1,S.*S.*n2);
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S^2 Q'_{ij}";
a=a+1;

[s_11, s_12, s_21, s_22] = exterior_product(n1, n2, n1_t, n2_t);
G(:,a) = 2*gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\partial_t Q_{ij}";
a=a+1;

s_11 = U.*Q_11_x + V.*Q_11_y;
s_12 = U.*Q_12_x + V.*Q_12_y; s_21 = s_12; s_22 = -s_11;
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_k \nabla_k Q_{ij}";
a=a+1;

s_11 = Omega_12.*Q_12;
s_12 =-Omega_12.*Q_11; 
s_21 =-Omega_12.*Q_11; 
s_22 =-Omega_12.*Q_12;
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\Omega_{ik} Q_{kj}";
a=a+1;


[s_11, s_12, s_21, s_22] = exterior_product(n1, n2, S.*U, S.*V);
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S u_i n_j";
a=a+1;

[s_11, s_12, s_21, s_22] = exterior_product(n1, n2, un.*S.*n1, un.*S.*n2);
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S u_k n_k Q'_{ij}";
a=a+1;

scalar = n1.*(n1_xx + n2_xy) + n2.*(n1_xy + n2_yy);
[s_11, s_12, s_21, s_22] = exterior_product(n1, n2, scalar.*n1, scalar.*n2);
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "(n_k \nabla_k S) Q'_{ij}";
a=a+1;

[s_11, s_12, s_21, s_22] = exterior_product(n1, n2, n1_xx+n2_xy, n1_xy+n2_yy );
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i \nabla_j S";
a=a+1;

[s_11, s_12, s_21, s_22] = exterior_product(S.*n1, S.*n2, B1, B2 );
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i S B_j";
a=a+1;

scalar = A_11.*(n1.^2 - n2.^2) + 2*A_12.*n1.*n2;
[s_11, s_12, s_21, s_22] = exterior_product( n1, n2, scalar.*n1, scalar.*n2 );
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "A'_{kl} Q'_{kl} Q'_{ij}";
a=a+1;

G(:,a) = gpu_integrate_symm_TF_tensor( A_11, A_12, A_12, -A_11, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "A'_{ij}";
a=a+1;

scalar = U.*B1 + V.*B2;
[s_11, s_12, s_21, s_22] = exterior_product( n1, n2, scalar.*n1, scalar.*n2 );
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_k B_k Q'_{ij}";
a=a+1;

%{
{["\eta = 0.022112, 0.24377 u_k \nabla_k Q_{ij}  +  -0.48715 S u_i n_j  +  0.47921 S u_k n_k Q'_{ij}  +  0.48773 u_k B_k Q'_{ij}  +  0.48553 u_i B_j  +  "                 ]}
[s_11, s_12, s_21, s_22] = exterior_product( U, V, B1, B2 );
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i B_j";
a=a+1;
%}

[s_11, s_12, s_21, s_22] = exterior_product( B1, B2, B1, B2 );
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B_i B_j";
a=a+1;

%{
(B_i B_j + B^2 Q_{ij} )' = 0

[s_11, s_12, s_21, s_22] = exterior_product( n1, n2, (B1.^2 + B2.^2).*n1, (B1.^2 + B2.^2).*n2 );
G(:,a) = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B^2 Q'_{ij}";
a=a+1;
%}



G_full = [G_full; G];
end %end weight iteration

G = G_full;
mask = mask0;
fprintf("Done integrating.\n")



function [s_11, s_12, s_21, s_22] = exterior_product(n1,n2,m1,m2)
  %symmetrized exterior product
  s_11 = n1.*m1;
  s_12 = (n1.*m2 + n2.*m1)/2;
  s_21 = s_12;
  s_22 = n2.*m2;
end


function vals = gpu_integrate_symm_TF_tensor( s_11, s_12, s_21, s_22, dx_vec, corners, size_vec, m, n, mask )
  %Only two degrees of freedom
  s_11 = (s_11 - s_22)/2;
  s_12 = (s_12 + s_21)/2;
  vals1 = gpu_integrate_data( s_11, [], dx_vec, corners, size_vec, m, n, mask );
  vals2 = gpu_integrate_data( s_12, [], dx_vec, corners, size_vec, m, n, mask );
  
  vals = [vals1; vals2];
  vals = gather(vals);
end