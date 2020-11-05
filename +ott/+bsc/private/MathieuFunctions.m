function F = Mathieu(x,r,q, type)

  if type == 'ce'
    F = Mathieu_ce(x, r, q);
  elseif type == 'se'
    F = Mathieu_se(x, r, q);
  else
    error('');
  end
end

function F = Mathieu_ce(x, r, q)
  %Angular Mathieu function type cosine with an approximations of N=100
  %terms, it has 3 arguments, the first one is the variable, the second
  %one is the eigenvalue and the third one is the ellipticity parameter.
  %Note: approximation can be improved by changing N.

  N=100;
  F=0;
  Ar = Fun_V_ce(r,q,N);


  if r==0 || mod(r,2)==0

    for s=1:N

      F = F + Ar(s)*cos(2*(s-1)*x);
    end

  else

    for s=1:N

      F = F + Ar(s)*cos((2*(s-1)+1)*x);
    end

  end
end


function F = Mathieu_se(x,r,q)
  %Angular Mathieu function type sin with an approximations of N=100
  %terms, it has 3 arguments, the first one is the variable, the second
  %one is the eigenvalue and the third one is the ellipticity parameter.
  %Note: approximation can be improved by changing N.

  N=100;
  F=0;
  B = Fun_V_se(r,q,N);
  F1=0;
  x0 = 0.1;

  if q==0

    F = sin(r*x);

  else

    if mod(r,2)==0

      for s=1:N

        F = F + B(s)*sin(2*s*x);
        F1 = F1 + B(s)*sin(2*s*x0);

      end

    else

      for s=1:N

        F = F + B(s)*sin((2*(s-1)+1)*x);
        F1 = F1 + B(s)*sin((2*(s-1)+1)*x0);

      end

    end

    if F1<0
      F=-F;
    end
  end
end

function M = Fun_MCnon(N,q)
  % Function that outputs the matrix M of odd coefficients for the calculation
  %of Mathieu functions type cosine. It carries 4 arguments; the first and the
  %second correspond to i and j respectively, the third argument is the
  %size of the matrix and the fourth argument is the parameter of ellipticity.

  for i=1:N    % rows
    for j=1:N % columns
      if j>=i
        if j==i
          if i==1
            M(i,j)=1+q;
          else
            M(i,j)=(2*(i-1) + 1)^2;
          end
        elseif j==i+1
          M(i,j)=q;
        else
          M(i,j)=0;
        end
      else
        M(i,j)=M(j,i);
      end
    end
  end
end


function M = Fun_MCpar(N,q)
  % Function that outputs the matrix M of even coefficients for the calculation
  %of Mathieu functions type cosine. It carries 4 arguments; the first and the
  %second correspond to i and j respectively, the third argument is the
  %size of the matrix and the fourth argument is the parameter of ellipticity.

  for i=1:N    % rows

    for j=1:N % columns

      if j>=i

        if j==i

          M(i,j)=(2*(i-1))^2;

        elseif j==i+1

          if i==1
            M(i,j)=sqrt(2)*q;
          else
            M(i,j)=q;
          end

        else
          M(i,j)=0;
        end

      else

        M(i,j)=M(j,i);
      end
    end
  end
end

function M = Fun_MSnon(N,q)
  % Function that outputs the matrix M of odd coefficients for the calculation
  %of Mathieu functions type sin. It carries 4 arguments; the first and the
  %second correspond to i and j respectively, the third argument is the
  %size of the matrix and the fourth argument is the parameter of ellipticity.

  for i=1:N    % rows
    for j=1:N % columns
      if j>=i
        if j==i
          if i==1
            M(i,j)=1-q;
          else
            M(i,j)=(2*(i-1) + 1)^2;
          end

        elseif j==i+1
          M(i,j)=q;

        else
          M(i,j)=0;
        end

      else
        M(i,j)=M(j,i);

      end
    end
  end
end


function M = Fun_MSpar(N,q)
  % Function that outputs the matrix M of even coefficients for the calculation
  %of Mathieu functions type sin. It carries 4 arguments; the first and the
  %second correspond to i and j respectively, the third argument is the
  %size of the matrix and the fourth argument is the parameter of ellipticity.

  for i=1:N    % rows

    for j=1:N % columns

      if j>=i

        if j==i

          M(i,j)=(2*i)^2;

        elseif j==i+1

          M(i,j)=q;

        else
          M(i,j)=0;
        end

      else

        M(i,j)=M(j,i);
      end
    end
  end
end

function  Ar = Fun_V_ce(r,q,N)
  % Calculate the k-th eigenvector of the angular Mathieu function
  % of type cosine with an approximation to order N

  if r==0
    M=Fun_MCpar(N,q);
    r1=1;

  elseif   mod(r,2)==0
    M=Fun_MCpar(N,q);
    r1=(r/2)+1;

  else
    M=Fun_MCnon(N,q);
    r1=(r+1)/2;
  end

  [v,d]=eig(M);   % Calculate eigenvectors, v, and eigenvalues d
  [~,j]=sort(diag(d));  % Sort by eigenvalues, j contains new positions
  A = v(:, j);

  %D=d(r1);    % eigenvalue a_r
  A=A(:,r1);    % eigenvector

  if r==0 || mod(r,2)==0
    A(1)= A(1)/sqrt(2);
  end

  Ar=A;

  if Ar(1)<0
    Ar=-Ar;
  end
end

function  B = Fun_V_se(r,q,N)
  % Calculate the k-th eigenvector of the angular Mathieu function
  % of type sin with an approximation to order N.

  if mod(r,2)==0

    M=Fun_MSpar(N,q);
    r1=r/2;

  else
    M=Fun_MSnon(N,q);
    r1=(r+1)/2;
  end

  [v,d]=eig(M);   % Calculate eigenvectors, v, and eigenvalues, d
  [~,j]=sort(diag(d));  % Sort by eigenvalues, j contains new positions
  B = v(:, j);

  %D=d(r1);    % eigenvalue a_r
  B=B(:,r1);    % eigenvector

  if B(1)<0
    B=-B;
  end
end
