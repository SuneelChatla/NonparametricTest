C   part of  VCsblcTest.f  by Stefan Sperlich created in 2018
C
C   Varying Coefficients Smoothed Backfitting Local Constant  Test
C
C   Test in varying coefficient models with cross sectional data 
C   loc.const. SB (has smoothing bias) for simulations 
C     for Mammen-Sperlich paper, submitted to Biometrika  
C
C   WITHOUT LINEAR PRIOR ESTIMATES FOR STARTING VALUES (so far)
C
C  this program can be used for H0: m_1 = constant
C     no matter whether m_1 is varying coeff or additive component (z_1=1) 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  SUBROUTINE  sbvctest  
C     loc-const smooth backfit estimation of varying coefficients mdodel
C  INPUTS :
C   x,z = designmatrices, 
C   y   = dependent observations
C   n   = number of observations
C   d   = number of dimensions 
C   h   = bandwidth vector
C   ng  = grid length (same for each dimension) 
C   maxit = maximal number of iterations (if maxi<1 then default=25)
C   conv= convergence cirteria: (if conv>0.1 then conv=0.01)
C           if  sum (xg_new - xg_old)^2 / sum xg_old^2 < conv STOP 
C   hd  = bandwidth for test statistic 
C   m   = number of trimming points for test on both sides of grid
C  OUTPUT :
C   to  	= test statistic
C   y     = (n) estimate of E[Y|X,Z] under H_0
C   maxit = number of iterations used
C   mf    = function estimates on grid under H_1
C   xg    = grid on which we estimated
C   m0    = estimated constant (under H_1)
      subroutine sbvctest(x,z,y,n,d,h,mf,xg,ng,maxit,conv,m,hd,to,m0)
	implicit real(kind=8) (a-h,o-z)
      integer n,d,it,i,j,k,l,ng,maxit,m
	real(kind=8) x(n,d),z(n,d),y(n),h(d),mie(ng,d,2),s(d),xg(ng,d)
     .            ,dens0(ng,d),dens(ng,ng,d*d),w(n,ng,d),dn
     .            ,m0,mf(ng,d),ymean,conv,m1(n),hd,to
      dn     = float(n)
      ymean  = sum(y)/dn
      if (maxit.LT.1) then 
        maxit = 25
      endif
      if (conv.GT.0.1) then 
        conv = 0.01
      endif
C      create grids for estimation, integration and interpolation
      do j=1,d
        s(j) = ( maxval(x(:,j))-minval(x(:,j)) ) / (ng-1.0)
        do i=1,ng,1	           
          xg(i,j) = minval(x(:,j)) + (i - 1.0)*s(j)     
        enddo
      enddo
C       calculate kernel weight matrices 
      do j=1,d,1
 	  w(:,:,j) = spread(x(:,j),2,ng)-spread(xg(:,j),1,n)
 	  w(:,:,j) = (w(:,:,j) / h(j)) ** 2.
        where (w(:,:,j).LT.1.0)  	  ! quartic kernel
          w(:,:,j) = 15./16. *((1.-w(:,:,j))**2.)/h(j)
        elsewhere
          w(:,:,j) = 0.0
        endwhere
        m1 =s(j)*(sum(w(:,2:(ng-1),j),2)+0.5*(w(:,1,j)+w(:,ng,j)))
        do k=1,ng	 ! projection kernel (must integrate to one)
          where (m1.GT.0.00001)  
            w(:,k,j) = w(:,k,j) / m1  * z(:,j)  ! times Z because VCM
          elsewhere
            w(:,k,j) = 0.0
          endwhere
        enddo
      enddo
C       now the marginal estimates for all dimensions
      mie = 0.0
      crit = 1.0
      i = 1
      do j=1,d,1
        dens0(:,j) = sum(w(:,:,j),1)
        m1(1:ng) = sum((w(:,:,j)*spread(z(:,j),2,ng)),1)
        where (m1(1:ng).GT.0.00001) 
          mie(:,j,1) = sum(w(:,:,j)*spread(y,2,ng),1)/m1(1:ng)  
          mie(:,j,2) = dens0(:,j) / m1(1:ng)
        endwhere
        do k=1,d,1
         if (j.NE.k) then
          do l=1,ng,1 
           where (m1(1:ng).GT.0.00001)
            dens(:,l,i)=sum(w(:,:,j)*spread(w(:,l,k),2,ng),1)/m1(1:ng)
           elsewhere
            dens(:,l,i) = 0.0
           endwhere
          enddo
         endif
         i = i+1
        enddo
      enddo	 
C      The backfitting iteration - first setting starting values
      w(1,:,:) = 1.0
      m0 = ymean
	mf = 0.0
	it = 0
      do while ((crit.GT.conv).AND.(maxit.GT.it)) 
	  it = it + 1
	  do j=1,d,1   ! loop over nonparametric functions
          mf(:,j) = mie(:,j,1) - m0 * mie(:,j,2)  
          do k=1,d,1
            if (k.NE.j) then   
              i = k+d*(j-1)
              m1(1:ng) = matmul(dens(:,2:(ng-1),i),mf(2:(ng-1),k))+
     *               .5*(mf(1,k)*dens(:,1,i)+mf(ng,k)*dens(:,ng,i))
              mf(:,j) = mf(:,j) - m1(1:ng)*s(k)  
            endif
          enddo
	  enddo
        m0 = ymean   ! the 2nd loop is for the constant
        do k=1,d,1
          m1(1) = sum( (dens0(2:(ng-1),k)*mf(2:(ng-1),k)) ) +
     *              .5*( mf(1,k)*dens0(1,k)+mf(ng,k)*dens0(ng,k) ) 
	    m0 = m0 - m1(1)*s(k) / n
        enddo
	  m1(1:d) = sum((mf-w(1,:,:))**2.,1)/(sum(w(1,:,:)**2.,1)+0.0001)
   	  crit = maxval(m1(1:d))
	  w(1,:,:) = mf
	enddo
      maxit = it
C      Test statistic for mf(:,1)	-- also with projection kernel 
 	w(:,:,1) = spread(x(:,1),2,ng)-spread(xg(:,1),1,n)
 	w(:,:,1) = (w(:,:,1) / hd ) ** 2.
      where (w(:,:,1).LT.1.0)  
        w(:,:,1) = 15./16. *((1.-w(:,:,1))**2.)/hd
      elsewhere
        w(:,:,1) = 0.0
      endwhere
      m1 =s(1)*(sum(w(:,2:(ng-1),1),2)+0.5*(w(:,1,1)+w(:,ng,1)))
      do k=1,ng	 ! now projection kernel (must integrate to one)
       where (m1.GT.0.00001)  
         w(:,k,1) = w(:,k,1) / m1   
       elsewhere
         w(:,k,1) = 0.0
       endwhere
      enddo
      dens0(:,2) = 0.0		! trimming in test statistic 
      dens0((1+m):(ng-m),2)=sum(w(:,(1+m):(ng-m),1),1)/dn    
!         we center the function before it enters the test statistic   
      ymean = ( sum((mf(2:(ng-1),1)*dens0(2:(ng-1),2)))+
     *          (mf(1,1)*dens0(1,2)+mf(ng,1)*dens0(ng,2))/2.)*s(1) 
      dens0(:,1) = mf(:,1) - ymean
      dens0(:,1) = dens0(:,1)**2.0 
      to = s(1)*( sum( (dens0(2:(ng-1),2)*dens0(2:(ng-1),1)) )+
     *          ( dens0(1,1)*dens0(1,2)+dens0(ng,1)*dens0(ng,2) )/2.)       
C     Interpolation on observed data for mf_2...mf_d ( mf_1=const under H_0 )
      y = ymean*z(:,1)   ! recall ymean=E[mf_1], needed for hatY under H_0
      mie(1:ng,2:d,1) = mf(:,2:d) 
      do j=2,d	! not for d=1 as this is in H_0 a constant	
        w(:,:,j) = spread(x(:,j),2,ng)-spread(xg(:,j),1,n)
        w(:,:,j) = w(:,:,j) / s(j)
        m1(1:(ng-1)) = mie(2:ng,j,1) - mie(1:(ng-1),j,1)        
        m1(ng)       = 0.0
        where ((w(:,:,j).GE.0.0).AND.(w(:,:,j).LT.1.0))
          w(:,:,j) =w(:,:,j)*spread(m1(1:ng),1,n)+spread(mie(:,j,1),1,n)
        elsewhere
          w(:,:,j) = 0.0
        endwhere
        w(:,1,j) = sum(w(:,:,j),2)
        where (x(:,j).LT.xg(1,j))   w(:,1,j) = mie(1,j,1)
        where (x(:,j).GT.xg(ng,j))  w(:,1,j) = mie(ng,j,1)
        y = y + w(:,1,j) *z(:,j)	! predictor	H_0, without constant m0
      enddo
      return
      end
CCCCCCCCCCCCCCCCCCCCCC END CCCCCCCCCCCCCCCCCCCCCCCCCC