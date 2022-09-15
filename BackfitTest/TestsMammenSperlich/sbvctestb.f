C   part of VCsblcTest.f  by Stefan Sperlich created in 2018
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
C  SUBROUTINE  sbvctestb   -- like vcs but with bootstrap  
C     loc-const smooth backfit estimation of varying coefficients mdodel
C  INPUTS :
C   x,z = designmatrices, 
C   y   = dependent observations
C   n   = number of observations
C   d   = number of dimensions 
C   h   = bandwidth vector
C   ng  = grid size (same for each dimension) 
C   maxit = maximum of iterations (if maxit=0 then default=25)
C   conv = convergence cirterion: (if conv>0.1 then conv=0.01)
C           if  sum (xg_new - xg_old)^2 / sum xg_old^2 < conv STOP  
C   hd  = bandwidth for test statistic 
C   rd  = random seed for replicability
C  OUTPUT :
C   pv  = bootstrap p-value
C
      subroutine sbvctestb(x,z,y,n,d,h,ng,mb,mit,cc,m,hd,to,nb,pv,rd) 
	implicit real(kind=8) (a-h,o-z)
      integer n,d,it,i,j,k,l,ng,mit,nb,m
	real(kind=8) x(n,d),z(n,d),y(n),h(d),mie(ng,d,2),s(d),xg(ng,d)
     .            ,ymean,dens0(ng,d),dens(ng,ng,d*d),w(n,ng,d),w1,w2 
     .            ,m0,mf(ng,d),cc,m1(ng),hd,to,pv,tb,m2(ng,d),pw,dn
     .	        ,rd,mb(n),densi(ng,2),yb(n)
      dn  = float(n) 
      if (mit.LT.1) then 
        mit = 25	 ! maxit
      endif
      if (cc.GT.0.1) then 
        cc = 0.01
      endif
C      create grids for estimation, integration and interpolation
      do j=1,d
        s(j) = ( maxval(x(:,j))-minval(x(:,j)) ) / (ng-1.0)
        do i=1,ng,1	           
          xg(i,j) = minval(x(:,j)) + (i - 1.0)*s(j)     
        enddo
      enddo
C      calculate density for test statistic, 
      w(:,:,1) = spread(x(:,1),2,ng)-spread(xg(:,1),1,n)
      w(:,:,1) = (w(:,:,1) / hd ) ** 2.
      where (w(:,:,1).LT.1.0)  
        w(:,:,1) = 15./16. *((1.-w(:,:,1))**2.)/hd
      elsewhere
        w(:,:,1) = 0.0
      endwhere
      yb =s(1)*(sum(w(:,2:(ng-1),1),2)+0.5*(w(:,1,1)+w(:,ng,1)))
      do k=1,ng	 ! now projection kernel (must integrate to one)
       where (yb.GT.0.00001)  
         w(:,k,1) = w(:,k,1) / yb   
       elsewhere
         w(:,k,1) = 0.0
       endwhere
      enddo
      densi(:,2) = 0.0		! trimming in test statistic 
      densi((1+m):(ng-m),2)=sum(w(:,(1+m):(ng-m),1),1)/dn    
C      now calculating the weights for the backfit
      do j=1,d,1
 	  w(:,:,j) = spread(x(:,j),2,ng)-spread(xg(:,j),1,n)
 	  w(:,:,j) = (w(:,:,j) / h(j)) ** 2.
        where (w(:,:,j).LT.1.0)  
          w(:,:,j) = 15./16. *((1.-w(:,:,j))**2.)/h(j)
        elsewhere
          w(:,:,j) = 0.0
        endwhere
        yb =s(j)*(sum(w(:,2:(ng-1),j),2)+0.5*(w(:,1,j)+w(:,ng,j)))
        do k=1,ng	 ! now projection kernel (must integrate to one)
          where (yb.GT.0.00001)  
            w(:,k,j) = w(:,k,j) / yb * z(:,j)
          elsewhere
            w(:,k,j) = 0.0
          endwhere
        enddo
      enddo
C       now the marginal estimates for all dimensions
      mie = 0.0
      i = 1
      do j=1,d,1
        dens0(:,j) = sum(w(:,:,j),1)
        m1 = sum((w(:,:,j)*spread(z(:,j),2,ng)),1)
        where (m1.GT.0.00001) 
          mie(:,j,2) = dens0(:,j) / m1
        endwhere
        do k=1,d,1
          if (j.NE.k) then
          do l=1,ng,1 
            where (m1.GT.0.00001)
              dens(:,l,i)=sum(w(:,:,j)*spread(w(:,l,k),2,ng),1)/m1
            elsewhere
              dens(:,l,i) = 0.0
            endwhere
          enddo
          endif
          i = i+1
        enddo
      enddo	 
C      golden cut Bootstrap
      pv = 0.0
      w1 = (1.-sqrt(5.))/2.
      w2 = (1.+sqrt(5.))/2.
      pw = (1.+sqrt(5.))/(2.*sqrt(5.))           
C       Bootstrap loop starts here
      mb = mb + sum(y-mb)/dn   ! center the prediction 
      S2P31M = 2147483647.0
      do lb=1,nb		! loop over bootstrap replicates
        SEED = rd	    ! generating uniform random vars
        DO 5 i=1,n
          SEED = MOD(16807.d0*SEED,S2P31M)
 5        yb(i) = SEED / 2.0**31.0
        rd = SEED 	 ! end of random generator 
   !     call ggubs(rd,n,yb)   
        where (yb.le.pw)
          yb = mb + (y-mb)*w1
        elsewhere
          yb = mb + (y-mb)*w2
        endwhere
      do j=1,d,1
        m1 = sum((w(:,:,j)*spread(z(:,j),2,ng)),1)
        where (m1.GT.0.00001) 
          mie(:,j,1) = sum(w(:,:,j)*spread(yb,2,ng),1) / m1 
        endwhere
      enddo	 
C      The backfitting iteration - first setting starting values
      ymean = sum(yb) / dn
      crit = 1.0
      it = 0
      m2 = 1.0
      m0 = ymean 
      mf = 0.0
      do while ((crit.GT.cc).AND.(mit.GT.it)) 
	  it = it + 1
	  do j=1,d,1
          mf(:,j) = mie(:,j,1) - m0 * mie(:,j,2)  
          do k=1,d,1
            if (k.NE.j) then   
              i = k+d*(j-1)
              m1 = matmul(dens(:,2:(ng-1),i),mf(2:(ng-1),k))+
     *               .5*(mf(1,k)*dens(:,1,i)+mf(ng,k)*dens(:,ng,i))
              mf(:,j) = mf(:,j) - m1*s(k)  
            endif
          enddo
	  enddo
        m0 = ymean
        do k=1,d,1
          m1(1) = sum( (dens0(2:(ng-1),k)*mf(2:(ng-1),k)) ) +
     *              .5*( mf(1,k)*dens0(1,k)+mf(ng,k)*dens0(ng,k) ) 
	    m0 = m0 - m1(1)*s(k) / n
        enddo
	  m1(1:d) = sum((mf-m2)**2.,1)/(sum(m2**2.,1)+0.0001)
   	  crit = maxval(m1(1:d))
	  m2 = mf
	enddo
C      Test statistic for mf(:,1)  center the function [ if required ] 
      densi(:,1) = mf(:,1) - (sum((mf(2:(ng-1),1)*densi(2:(ng-1),2)))+
     *          (mf(1,1)*densi(1,2)+mf(ng,1)*densi(ng,2))/2.)*s(1) ! center
      densi(:,1) = densi(:,1)**2.0 
      tb = s(1)*(sum( (densi(2:(ng-1),2)*densi(2:(ng-1),1)) )+
     *          ( densi(1,1)*densi(1,2)+densi(ng,1)*densi(ng,2))/2.) 
      if (tb.GT.to)  pv = pv+1.0
      enddo ! bootstrap loop
      pv = pv / float(nb)	 ! p-value
      return
      end
CCCCCCCCCCCCCCCCCCCCCC END CCCCCCCCCCCCCCCCCCCCCCCCCC