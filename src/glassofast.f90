      subroutine glassofast(n,S,L,thr,maxIt,msg,warm,X,W,info)
!
!     .. Scalar Arguments ..
      implicit double precision(a-h, o-z)
      integer n, warm, msg, maxIt, iter
!     ..
!     .. Array Arguments ..
      double precision S(n,n), L(n,n), X(n,n), W(n,n)
!     ..
!
!  Purpose
!  =======
!
!  This subroutine computes the L1 regularized covariance matrix estimate
!  using the algorithm described in the paper:
!    J. Friedman, T. Hastie, R. Tibshirani:
!    Sparse inverse covariance estimation with the graphical lasso
!    Biostatistics, 9(3):432-441, July 2008.
!  This implementation is documented in the paper:
!    M. A. Sustik B. Calderhead:
!    Efficient Glasso implementation
!    Journal of Machine Learning Research, submitted.
!
!  Arguments
!  =========
!
!  n      (input) integer
!         The dimension of the input matrix s
!
!  S      (input) double precision array, dimension n x n
!         The empirical covariance matrix
!
!  L      (input) double precision array, dimension n x n
!         Regularization matrix (symmetric)
!
!  thr    (input) double precision
!         Convergence threshold
!
!  maxIt  (input) integer
!         Maximum number of whole matrix sweeps
!
!  msg    (input) integer
!         Controls amount of messages printed
!
!  warm   (input) integer
!         flag indicating cold versus warm start, see also X, W
!
!  X      (input/output) double precision array, dimension n x n
!         Inverse covariance matrix estimate
!
!  W      (input/output) double precision array, dimension n x n
!         Covariance matrix estimate
!
!  info   (output) integer
!         Indicates errors

double precision, dimension (:), allocatable :: Wd, WXDj
double precision EPS
parameter (EPS = 1.1e-16)
allocate(Wd(1:n), stat = info)
allocate(WXDj(1:n), stat = ierr)
info = info + ierr
if (info .ne. 0) then
   return
endif
shr = sum(abs(S))
do i = 1,n
   shr = shr - abs(S(i, i))
enddo
if (shr .eq. 0.0) then
!  S is diagonal.
   W = 0.0
   X = 0.0
   do i = 1,n
      W(i,i) = W(i,i) + L(i,i)
   enddo
   X = 0.0
   do i = 1,n
      X(i,i) = 1.0/max(W(i,i),eps)
   enddo
   return
endif
shr = thr*shr/(n-1)
if (warm .eq. 0) then
   W = S
   X = 0.0
   do i = 1,n
      X(i,i) = 1.0
   enddo
else
   do i = 1,n
     tmp = X(i,i)
     X(1:n,i) = X(1:n,i)/tmp
  end do
endif
do i = 1,n
   W(i,i) = S(i,i) + L(i,i)
   Wd(i) = W(i,i)
enddo
do iter = 1,maxIt
!   JC: print iterations to the console (ok with CRAN)
    if (msg .ne. 0)  call intpr('iter:',-1,iter,1)
dw = 0.0
   do j = 1,n
      WXDj(1:n) = 0.0
!     We exploit sparsity of X when computing column j of W*X*D:
      do i = 1,n
         if (X(i,j) .ne. 0.0) then
            WXDj = WXDj - W(:,i)*X(i,j)
         endif
      enddo
      do
         dlx = 0.0
         do i = 1,n
            if (i .ne. j) then
               diff = S(i,j) - WXDj(i) - W(i,j) - Wd(i)*X(i,j)
               tmp = abs(diff) - L(i,j)
               u = -X(i,j)
               X(i,j) = 0.0
               if (tmp .gt. 0.0) then
                  X(i,j) = -sign(tmp,diff)/Wd(i)
                  u = u + X(i,j)
               endif
               if (u .ne. 0.0) then
                  WXDj(1:n) = WXDj(1:n) - W(:,i)*u
                  dlx = max(dlx, abs(u))
               endif
            endif
         enddo
         if (dlx .lt. thr) then
            exit
         endif
      enddo
      dw = max(dw,sum(abs(WXDj(1:n))) - abs(WXDj(j)))
      WXDj(1:n) = WXDj(1:n) + W(:,j)
      W(:,j) = WXDj(1:n)
      W(j,:) = WXDj(1:n)
      W(j,j) = Wd(j)
   enddo
   if (dw .le. shr) then
      exit
   endif
enddo
do i = 1,n
   tmp = 1/(sum(X(:,i)*W(:,i)))
   X(1:n,i) = tmp*X(1:n,i)
enddo
do i = 1,n
   X(1:n,i) = (X(1:n,i) + X(i,1:n))/2;
   X(i, 1:n) = X(1:n,i) 
enddo
maxIt = iter ! JC : Update the number of iterations
return
end subroutine glassofast
