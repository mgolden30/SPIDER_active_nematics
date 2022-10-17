addpath('../bin')
fprintf('Starting vector regression without curl... (momentum library)\n');

num_lib     = 55;
num_windows = 10*num_lib;
G = zeros(2*num_windows, num_lib);

labels = cell(num_lib, 1);
%Loop over windows and integrate.
%{
  We will be integrating vector fields to get G using
  G = \int dV \int dt grad(phi) cross (Library vector)
  We use the 2D cross product: (a,b)cross(c,d) = ad-bc.
  This integration removes dependence on pressure, since the integral will
  vanish identically.
%}

m = 4*[1 1 1];
n = m;

corners = zeros(3, num_windows);

size_vec= round( 2./dx_vec); %integrate over a few characteristic length and time scales

edge = 12;
for i=1:num_windows
    corners(:,i) = [ randi(Ly-size_vec(1)-2*edge) randi(Lx-size_vec(2)-2*edge) randi(Lt-size_vec(3)-2*edge) ] + [1 1 1]*edge;
end

%Here are the algebraically independent derivatives
un = U.*n1 + V.*n2;

D = U_x + V_y;
A_11 =  (U_x - V_y)/2;
A_22 = -A_11;
A_12 =  (U_y + V_x)/2;
A_21 =  A_12;

Ann = A_11.*( n1.^2 - n2.^2 ) + 2*A_12.*n1.*n2;

S  = n1_x + n2_y;
B1 = -(n1.*n1_x + n2.*n1_y);
B2 = -(n1.*n2_x + n2.*n2_y);
Omega_12 = (U_y - V_x)/2;

S_x = n1_xx + n2_xy; S_y = n1_xy + n2_yy;
[B1_x, B1_y, ~] = gradient3d(B1, dx_vec);
[B2_x, B2_y, ~] = gradient3d(B2, dx_vec);

[U_yy, U_xx, U_tt, U_xy, U_xt, U_yt] = secondgradient3d( U, dx_vec );
[V_yy, V_xx, V_tt, V_xy, V_xt, V_yt] = secondgradient3d( V, dx_vec );





G_full = [];
lengths = size_vec.*dx_vec;
weights = { 0*x+1, cos( pi*x/lengths(2) ), cos( pi*y/lengths(1)), cos( pi*t/lengths(3)) };
mask0 = mask;

for ii = 1:numel(weights)
a = 1;

clear functions            %We will recompute derivatives of the mask
mask = mask0.*weights{ii}; %Multiply your weight by anything you want




G(:,a) = gpu_integrate_vector_curl( S.*n1, S.*n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S n_i";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( U, V, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( un.*n1, un.*n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i n_j u_j";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( U_t, V_t, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\partial_t u_i";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( n1.*(n1.*U_t + n2.*V_t), n2.*(n1.*U_t + n2.*V_t), dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i n_j \partial_t u_j";
a = a+1;

%{
remove since div(Q) + B - Sn = 0.
G(:,a) = gpu_integrate_vector_curl( B1, B2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B_i";
a = a+1;
%}

G(:,a) = gpu_integrate_vector_curl( U.*S.*S, V.*S.*S, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i S^2";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( n1.*un.*S.*S, n2.*un.*S.*S, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i n_j u_j S^2";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( un.*n1_t, un.*n2_t, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_j n_j \partial_t n_i";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( n1.*(U.*n1_t + V.*n2_t), n2.*(U.*n1_t + V.*n2_t), dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i u_j \partial_t n_j";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( S.*un.*B1, S.*un.*B2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S n_j u_j B_i";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( S.*(U.*B1 + V.*B2).*n1, S.*(U.*B1 + V.*B2).*n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S n_i u_j B_j";
a = a+1;

% S \nabla u

G(:,a) = gpu_integrate_vector_curl( S.*Ann.*n1, S.*Ann.*n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S A'_{jk} Q'_{jk} n_i";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( S.*(n1.*U_x + n2.*V_x), S.*(n1.*U_y + n2.*V_y), dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S n_j \nabla_i u_j";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( S.*(n1.*U_x + n2.*U_y), S.*(n1.*V_x + n2.*V_y), dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S n_j \nabla_j u_i";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( n1.*S.*(U_x + V_y), n2.*S.*(U_x + V_y), dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S n_j \nabla_j u_i";
a = a+1;

% u \nabla S

G(:,a) = gpu_integrate_vector_curl( n1.*un.*(n1.*S_x + n2.*S_y), n2.*un.*(n1.*S_x + n2.*S_y), dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i (n_j u_j)( n_k \nabla_k S)";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( U.*(n1.*S_x + n2.*S_y), V.*(n1.*S_x + n2.*S_y), dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i n_j \nabla_j S";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( un.*S_x, un.*S_y, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_j n_j \nabla_i S";
a = a+1; 

% \nabla \nabla u

v1 = n1.*n1.*U_xx + n2.*n2.*U_yy + 2*n1.*n2.*(U_xy);
v2 = n1.*n1.*V_xx + n2.*n2.*V_yy + 2*n1.*n2.*(V_xy);
G(:,a) = gpu_integrate_vector_curl( v1, v2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_j n_k \nabla_j \nabla_k u_i";
a = a+1;

v1 = n1.*n1.*U_xx + n2.*n2.*V_xy + n1.*n2.*(U_xy + V_xx);
v2 = n1.*n1.*U_xy + n2.*n2.*V_yy + n1.*n2.*(U_yy + V_xy);
G(:,a) = gpu_integrate_vector_curl( v1, v2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_j n_k \nabla_i \nabla_j u_k";
a = a+1;

% u B B

G(:,a) = gpu_integrate_vector_curl( U.*(B1.^2 + B2.^2), V.*(B1.^2 + B2.^2), dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B^2 u_i";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( un.*(B1.^2 + B2.^2).*n1, un.*(B1.^2 + B2.^2).*n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B^2 u_j n_j n_i";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( B1.*(U.*B1 + V.*B2), B2.*(B1.*U + B2.*V).*n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B_i B_j u_j";
a = a+1;

%B \nabla u

G(:,a) = gpu_integrate_vector_curl( B1.*Ann, B2.*Ann.*n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B_i A'_{jk} Q'_{jk}";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( B1.*(U_x + V_y), B2.*(U_x + V_y), dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B_i \nabla_j u_k";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( B1.*U_x + B2.*U_y, B1.*V_x + B2.*V_y, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B_j \nabla_j u_i";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( B1.*U_x + B2.*V_x, B1.*U_y + B2.*V_y, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B_j \nabla_i u_j";
a = a+1;

% u \nabla B

G(:,a) = gpu_integrate_vector_curl( un.*(n1.*B1_x + n2.*B1_y), un.*(n1.*B2_x + n2.*B2_y), dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_j u_j n_k \nabla_k B_i";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( U.*(B1_x + B2_y), V.*(B1_x + B2_y), dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i \nabla_j B_j";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( U.*B1_x + V.*B2_x, U.*B1_y + V.*B2_y, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_j \nabla_i B_j";
a = a+1;

G(:,a) = gpu_integrate_vector_curl( U.*B1_x + V.*B1_y, U.*B2_x + V.*B2_y, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_j \nabla_j B_i";
a = a+1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stress Tensors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s_xx = n1.*n1; s_xy = n1.*n2; s_yx = n2.*n1; s_yy = n2.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i n_j";
a=a+1;

s_xx = S.*S.*n1.*n1; s_xy = S.*S.*n1.*n2; s_yx = s_xy; s_yy = S.*S.*n2.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S^2 n_i n_j";
a=a+1;

%S u

s_xx = S.*n1.*U; s_xy = S.*n1.*V; s_yx = S.*n2.*U; s_yy = S.*n2.*V;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S n_i u_j";
a=a+1;

s_xx = S.*n1.*U; s_xy = S.*n2.*U; s_yx = S.*n1.*V; s_yy = S.*n2.*V;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S n_j u_i";
a=a+1;

scalar = S.*un;
s_xx = scalar.*n1.*n1; s_xy = scalar.*n1.*n2; s_yx = s_xy; s_yy = scalar.*n2.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S u_k n_k Q'_{ij}";
a=a+1;

% \nabla S
scalar = n1.*S_x + n2.*S_y;
s_xx = scalar.*n1.*n1; s_xy = scalar.*n1.*n2; s_yx = s_xy; s_yy = scalar.*n2.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "(n_k \nabla_k S) Q'_{ij}";
a=a+1;

s_xx = n1.*S_x; s_xy = n1.*S_y; s_yx = n2.*S_x; s_yy = S_y.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i \nabla_j S";
a=a+1;

s_xx = n1.*S_x; s_xy = n2.*S_x; s_yx = n1.*S_y; s_yy = S_y.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_j \nabla_i S";
a=a+1;

% S B

s_xx = n1.*S.*B1; s_xy = n1.*S.*B2; s_yx = n2.*S.*B1; s_yy = S.*n2.*B2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i S B_j";
a=a+1;

s_xx = n1.*S.*B1; s_xy = n2.*S.*B1; s_yx = n1.*S.*B2; s_yy = S.*n2.*B2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_j S B_i";
a=a+1;

% \nabla u
scalar = U_x + V_y;
s_xx = scalar.*n1.*n1; s_xy = scalar.*n1.*n2; s_yx = s_xy; s_yy = scalar.*n2.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "(\nabla_k u_k) Q'_{ij}";
a=a+1;

scalar = Ann;
s_xx = scalar.*n1.*n1; s_xy = scalar.*n1.*n2; s_yx = s_xy; s_yy = scalar.*n2.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "A'_{kl} Q'_{kl} Q'_{ij}";
a=a+1;

s_xx = (n1.*U_x + n2.*V_x).*n1; s_xy = (n1.*U_x + n2.*V_x).*n2; s_yx = (n1.*U_y + n2.*V_y).*n1; s_yy = (n1.*U_y + n2.*V_y).*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "(\nabla_i u_k) n_k n_j";
a=a+1;

s_xx = (n1.*U_x + n2.*V_x).*n1; s_xy = (n1.*U_y + n2.*V_y).*n1; s_yx = (n1.*U_x + n2.*V_x).*n2; s_yy = (n1.*U_y + n2.*V_y).*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "(\nabla_j u_k) n_k n_i";
a=a+1;

s_xx = A_11; s_xy = A_12; s_yx = A_12; s_yy = A_22;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "A'_{ij}";
a=a+1;

%{
A_{ij} - \Omega_{ij} is curlless

s_xx = 0*A_11; s_xy = Omega_12; s_yx = -Omega_12; s_yy = 0*A_11;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\Omega_{ij}";
a=a+1;
%}

% u B

scalar = U.*B1 + V.*B2;
s_xx = scalar.*n1.*n1; s_xy = scalar.*n1.*n2; s_yx = s_xy; s_yy = scalar.*n2.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_k B_k Q'_{ij}";
a=a+1;
  
s_xx = un.*B1.*n1; s_xy = un.*B1.*n2; s_yx = un.*B2.*n1; s_yy = un.*B2.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_k n_k B_i n_j";
a=a+1;

s_xx = un.*B1.*n1; s_xy = un.*B2.*n1; s_yx = un.*B1.*n2; s_yy = un.*B2.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_k n_k B_j n_i";
a=a+1;

%{
%Remove these two library terms with identity
%{["\eta = 4.1398e-15, -0.57735 u_k B_k Q'_{ij}  +   0.57735 u_k n_k B_i n_j  +  -0.57735 u_j B_i  +  "     ]}
%{["\eta = 5.6537e-15,  0.57735 u_k B_k Q'_{ij}  +  -0.57735 u_k n_k B_j n_i  +   0.57735 u_i B_j  +  "      ]}


s_xx = U.*B1; s_xy = U.*B2; s_yx = V.*B1; s_yy = V.*B2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i B_j";
a=a+1;

s_xx = U.*B1; s_xy = V.*B1; s_yx = U.*B2; s_yy = V.*B2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_j B_i";
a=a+1;
%}

% \nabla B

scalar = B1_x + B2_y;
s_xx = scalar.*n1.*n1; s_xy = scalar.*n1.*n2; s_yx = s_xy; s_yy = scalar.*n2.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "(\nabla_k B_k) Q'_{ij}";
a=a+1;

s_xx = (n1.*B1_x + n2.*B1_y).*n1; 
s_xy = (n1.*B1_x + n2.*B1_y).*n2; 
s_yx = (n1.*B2_x + n2.*B2_y).*n1; 
s_yy = (n1.*B2_x + n2.*B2_y).*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "(n_k \nabla_k B_i) n_j";
a=a+1;

s_xx = (n1.*B1_x + n2.*B1_y).*n1; 
s_xy = (n1.*B2_x + n2.*B2_y).*n1; 
s_yx = (n1.*B1_x + n2.*B1_y).*n2; 
s_yy = (n1.*B2_x + n2.*B2_y).*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "(n_k \nabla_k B_j) n_i";
a=a+1;


s_xx = B1_x; 
s_xy = B2_x; 
s_yx = B1_y; 
s_yy = B2_y;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\nabla_i B_j";
a=a+1;

s_xx = B1_x; 
s_xy = B1_y; 
s_yx = B2_x; 
s_yy = B2_y;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\nabla_j B_i";
a=a+1;


%From BB

scalar = B1.^2 + B2.^2;
s_xx = scalar.*n1.*n1; s_xy = scalar.*n1.*n2; s_yx = s_xy; s_yy = scalar.*n2.*n2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B^2 Q'_{ij}";
a=a+1;

%{
s_xx = B1.*B1; s_xy = B1.*B2; s_yx = s_xy; s_yy = B2.*B2;
G(:,a) = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B_i B_j";
a=a+1;
%}


G_full = [G_full; G];
end %end weight iteration

G = G_full;
mask = mask0;
fprintf("done integrating.\n")









function vals = gpu_integrate_vector_curl( U, V, dx_vec, corners, size_vec, m, n, mask )
  %Integrates U_y - V_x with integration by parts.
  vals1 =        gpu_integrate_data( U, [], dx_vec, corners, size_vec, m, n, mask );
  vals2 =        gpu_integrate_data( V, [], dx_vec, corners, size_vec, m, n, mask );
  %  vals = vals - gpu_integrate_data( V, [2], dx_vec, corners, size_vec, m, n, mask );
  vals = [vals1; vals2];
end

function vals = gpu_integrate_stress_tensor( s_xx, s_xy, s_yx, s_yy, dx_vec, corners, size_vec, m, n, mask )
  %integrates the vector field \nabla_j\sigma_{ij} in the manner described
  %in the previous function
  vals1 =         gpu_integrate_data( s_xx, [2], dx_vec, corners, size_vec, m, n, mask );
  vals1 =  vals1+ gpu_integrate_data( s_xy, [1], dx_vec, corners, size_vec, m, n, mask );
  vals2 =         gpu_integrate_data( s_yx, [2], dx_vec, corners, size_vec, m, n, mask );
  vals2 =  vals1+ gpu_integrate_data( s_yy, [1], dx_vec, corners, size_vec, m, n, mask );
  vals = [vals1;vals2];
%  vals = vals + gpu_integrate_data( s_xy,        [1,1], dx_vec, corners, size_vec, m, n, mask );
%  vals = vals - gpu_integrate_data( s_yx,        [2,2], dx_vec, corners, size_vec, m, n, mask );
end