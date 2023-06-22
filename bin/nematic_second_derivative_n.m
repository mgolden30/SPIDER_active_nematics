function [n1_xx, n1_xy, n1_yy, n2_xx, n2_xy, n2_yy] = nematic_second_derivative_n( n1, n2, dx_vec )
  n1_xx    = 0*n1;
  n1_xy    = 0*n1;
  n1_yy    = 0*n1;
  n2_xx    = 0*n1;
  n2_xy    = 0*n1;
  n2_yy    = 0*n1;
  
  %sum      = 0*nx_x + 0*ny_y;
  %eighth_term_1 = 0*sum;
  %eighth_term_2 = 0*sum;

  %https://www.mathworks.com/matlabcentral/answers/16996-assign-multiple-variables  
  [dy dx dt] = feval(@(x) x{:}, num2cell(dx_vec));

  dim = size(n1);
  for j = 2:dim(1)-1
    for i = 2:dim(2)-1
      for k = 2:dim(3)-1
        %copy a stencil from nx and ny
        center = [n1(j,i,k) n2(j,i,k)];
        top    = [n1(j,i+1,k) n2(j,i+1,k)];
        bottom = [n1(j,i-1,k) n2(j,i-1,k)];
        left   = [n1(j-1,i,k) n2(j-1,i,k)];
        right  = [n1(j+1,i,k) n2(j+1,i,k)];
  
        tr = [n1(j+1,i+1,k) n2(j+1,i+1,k)];
        tl = [n1(j-1,i+1,k) n2(j-1,i+1,k)];
        br = [n1(j+1,i-1,k) n2(j+1,i-1,k)];
        bl = [n1(j-1,i-1,k) n2(j-1,i-1,k)];

        %Now we will flip vectors as needed.
        %center*top' is a dot product.
        if center*top' < 0
          top = -top;
        end
        if center*bottom' < 0
          bottom = -bottom;
        end
        if center*right' < 0
          right = -right;
        end
        if center*left' < 0
          left = -left;
        end
        if center*tr' < 0
          tr = -tr;
        end
        if center*tl' < 0
          tl = -tl;
        end
        if center*br' < 0
          br = -br;
        end
        if center*bl' < 0
          bl = -bl;
        end
        
        n1_xx(j,i,k) = ( bottom(1) - 2*center(1) + top(1)   )/(dx*dx);
        n1_xy(j,i,k) = ( tr(1) - tl(1) - br(1) + bl(1)      )/(4*dx*dy);
        n1_yy(j,i,k) = ( left(1)   - 2*center(1) + right(1) )/(dy*dy);
        
        n2_xx(j,i,k) = ( bottom(2) - 2*center(2) + top(2)   )/(dx*dx);
        n2_xy(j,i,k) = ( tr(2) - tl(2) - br(2) + bl(2)      )/(4*dx*dy);
        n2_yy(j,i,k) = ( left(2)   - 2*center(2) + right(2) )/(dy*dy);
      end %k
    end %i
  end %j
end %func
