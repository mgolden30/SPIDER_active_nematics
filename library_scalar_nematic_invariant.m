addpath('../bin')

fprintf('Starting scalar regression...\n');

num_lib         = 9;
num_windows     = num_lib*10;
nw = num_windows; %some shorthand
nl = num_lib;

labels = cell(nl,1);
corners = zeros(3, num_windows);
size_vec= round( 2./dx_vec );

%size_vec = [72    72    50]
Ly = size(U,1);
Lx = size(U,2);
Lt = size(U,3);
for i=1:num_windows
  L = 10;
  corners(:,i) = [ randi(Ly-size_vec(1)-2*L) randi(Lx-size_vec(2)-2*L) randi(Lt-size_vec(3)-2*L) ] + [L L L];
end

%The polynomial weights uses these powers: x^m(1-x)^n
m = 4*[1 1 1];
n = m;

G = zeros(num_windows, num_lib);

%Helpful shorthands
un = U.*n1 + V.*n2;
B1 = -(n1.*n1_x + n2.*n1_y);
B2 = -(n1.*n2_x + n2.*n2_y);
S  = n1_x + n2_y;



G_full = [];
lengths = size_vec.*dx_vec;
%weights = { 0*x+1, cos( pi*x/lengths(2) ), cos( pi*y/lengths(1)), cos( pi*t/lengths(3)) };
mask0 = mask;

A_22 = -A_11;

for ii = 1:1%numel(weights)
a = 1;

clear functions            %We will recompute derivatives of the mask
mask = mask0;%.*weights{ii}; %Multiply your weight by anything you want






G(:,a) = gpu_integrate_data( 0*U+1, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "1";
a = a+1;

G(:,a) = gpu_integrate_data( S.^2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S^2";
a = a+1;

G(:,a) = gpu_integrate_data( S.*un, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i n_i S";
a = a+1;

G(:,a) =          gpu_integrate_data( S.*n1, [2], dx_vec, corners, size_vec, m, n, mask );
G(:,a) = G(:,a) + gpu_integrate_data( S.*n2, [1], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\nabla_i (S n_i)";
a = a+1;

G(:,a) =          gpu_integrate_data( U, [2], dx_vec, corners, size_vec, m, n, mask );
G(:,a) = G(:,a) + gpu_integrate_data( V, [1], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\nabla_i u_i";
a = a+1;

Ann = A_11.*n1.*n1 + 2*A_12.*n1.*n2 + A_22.*n2.*n2;
G(:,a) = gpu_integrate_data( Ann, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "A'_{ij} n_i n_j";
a = a+1;

G(:,a) = gpu_integrate_data( U.*B1 + V.*B2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i B_i";
a = a+1;

G(:,a) =          gpu_integrate_data( B1, [2], dx_vec, corners, size_vec, m, n, mask );
G(:,a) = G(:,a) + gpu_integrate_data( B2, [1], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\nabla_i B_i";
a = a+1;

B1 = n1.*n1_x + n2.*n1_y;
B2 = n1.*n2_x + n2.*n2_y;
G(:,a) = gpu_integrate_data( B1.^2 + B2.^2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "B^2";
a = a+1;


G(:,a) = gpu_integrate_data( sqrt( U_x.^2 + U_y.^2 + V_x.^2 + V_y.^2 ), [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\sqrt{(\nabla_i u_j)(\nabla_i u_j)}";
a = a+1;

G_full = [G_full; G];
end %end weight iteration

G = G_full;
mask = mask0;
fprintf("done integrating.\n")