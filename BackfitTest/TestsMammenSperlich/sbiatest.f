C  Stefan Sperlich, created in 2018  (sobroutine backfon from mamt2.f)
C
C  this file contains only the subroutine for the main estimation
C
C     All X are on similar scale around interval [-ra/2;ra/2]
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      SUBROUTINE  sbiatest	for loc.lin. smooth backfit
C                  with interaction m_12 and its decomposition
C                  Calculates test statistic value
C  INPUTS :
C   x   = designmatrix, all covariates rescaled to about [-ra/2;ra/2] 
C   ra  = range length  
C   y   = dependent observations	-- will be changed, see below !
C   h   = bandwidth vector for one-dim components
C   hd  = bandwidth vector for interaction m_12 (called e_12)
C   n   = number of observations
C   d   = number of dimensions
C   ng  = grid size in each direction
C   xg  = grid on X on which integral is calculated & functions estimated
C   trim = threshold for trimming (see papers with Dette and Tjostheim)
C  OUTPUT :
C   y  = regression estimates under H_0 without constant
C   to = test statistic value
C   m3 = function estimates (for x3,... without x1,x2)  
C   m12= function estimates (for x1,x2 under H_0 decomposition) 
      subroutine sbiatest(x,n,d,h,hd,ng,y,to,trim,ra,m3,m12)  
      implicit real(kind=8) (a-h,o-z)
      integer  n,d,it,i,j,k,l,ng,maxit,jk,jl,gn,dd 
      parameter(maxit=15,concr=0.00001,tol=0.0001)
      real(kind=8)  x(n,d),y(n),xg(ng),h(d),hd(2),s,ra  
     .    ,d11(ng*ng,ng,d-2),d21(ng,ng,ng,d-2),d31(ng,ng,d-2,d-3)
     .    ,d12(ng*ng,ng,d-2),d22(ng,ng,ng,d-2),d32(ng,ng,d-2,d-3)
     .    ,d13(ng*ng,ng,d-2),d23(ng,ng,ng,d-2),d33(ng,ng,d-2,d-3)
     .    ,d14(ng*ng,ng,d-2),d24(ng,ng,ng,d-2),d34(ng,ng,d-2,d-3)
     .    ,d15(ng*ng,ng,d-2),d25(ng,ng,ng,d-2)
     .    ,d16(ng*ng,ng,d-2),d26(ng,ng,ng,d-2) 
     .    ,co(ng*ng,6),df(ng,2,d-2),m0(ng),m1(ng),m2(ng),wo(n,ng,d)
     .    ,m12(ng*ng),md12(ng*ng,2),m3(ng,d-2),md3(ng,d-2),to,trim     
      real*8, ALLOCATABLE:: w(:,:,:),mde(:,:),ed3(:,:),e12(:),e3(:,:)
 ! 777  format(1X,3F10.3) 
        gn = ng*ng
        dd = d-2
      dn = float(n)
      cy = sum(y,1) / dn
      y = y - cy			! look at centered responses
C       creation of grid for estimation and integration
      s = ra / (ng-1.0)    
      do j=1,ng,1	           
        xg(j) = (j - 1.0)*s  - 0.5*ra  
      enddo
C        calculating the kernel weights
      ALLOCATE (w(n,ng,d),e12(n))
      do j=1,d,1	           
        wo(:,:,j) = spread(x(:,j),2,ng)-spread(xg(:),1,n)
        w(:,:,j) = (wo(:,:,j) / h(j)) ** 2.
        where (w(:,:,j).LT.1.0)  
          w(:,:,j) = 0.9375*((1.-w(:,:,j))**2.) / h(j)
        elsewhere
          w(:,:,j) = 0.0
        endwhere	  ! now for projection kernel
        e12 = s*(sum(w(:,2:(ng-1),j),2)+0.5*(w(:,1,j)+w(:,ng,j)))
        do k=1,ng
         where (e12.GT.tol)  
           w(:,k,j) = w(:,k,j) / e12 
         elsewhere
           w(:,k,j) = 0.0
         endwhere
        enddo
      enddo
      DEALLOCATE (e12)
C      creation of auxiliary terms starting with nuissance directions 
      do j=3,d
       it = j-2
       m1 = sum(w(:,:,j),1)			 !  for p_3
       m2 = sum(w(:,:,j)*wo(:,:,j),1)	 !  for p_3^3  hilfsf
       where (m1.GT.tol)
         m3(:,it)  = sum(w(:,:,j)*spread(y,2,ng),1) / m1
         df(:,1,it) = m2 / m1
       elsewhere
         m3(:,it)   = 0.0	 !	m_3
         df(:,1,it) = 0.0	 !
       endwhere
       where (abs(m2).GT.tol)
         md3(:,it) = sum(w(:,:,j)*wo(:,:,j)*spread(y,2,ng),1) / m2
         df(:,2,it) = sum(w(:,:,j)*wo(:,:,j)**2.,1) / m2
       elsewhere
         md3(:,it)  = 0.0	 !	m_3^3  hilfsf
         df(:,2,it) = 0.0   !    
       endwhere
       do i=1,ng 
        if (m1(i).GT.tol) then
         l=1
         do k=3,d
          if (j.NE.k) then
           d31(i,:,it,l) = sum(w(:,:,k)*
     *                           spread(w(:,i,j),2,ng),1) / m1(i)
           d32(i,:,it,l) = sum(w(:,:,k)*wo(:,:,k)*
     *                           spread(w(:,i,j),2,ng),1) / m1(i)
           l = l+1
          endif
         enddo ! k
        endif
        if (abs(m2(i)).GT.tol) then
         l=1
         do k=3,d
          if (j.NE.k) then
           d33(i,:,it,l) = sum(w(:,:,k)*
     *                  spread(w(:,i,j)*wo(:,i,j),2,ng),1) / m2(i)
           d34(i,:,it,l) =  sum(w(:,:,k)*wo(:,:,k)*
     *                  spread(w(:,i,j)*wo(:,i,j),2,ng),1) / m2(i)
           l = l+1
          endif
         enddo ! k
        endif
       enddo ! i 
       d23(:,:,:,it) = spread(spread(m1,1,ng),2,ng) ! p_3
       d26(:,:,:,it) = spread(spread(m2,1,ng),2,ng) ! p_3^3
      enddo ! j
      co = 0.0
      d1 = 0.0     ! now for the interaction of x1, x2 
      do i=1,ng,1
        k = 1+(i-1)*ng
        l = i*ng
        m1 = sum(w(:,:,1)*spread(w(:,i,2),2,ng),1)  ! p_12
        m2 = sum(w(:,:,1)*wo(:,:,1)*spread(w(:,i,2),2,ng),1) ! p_12^1
        m0 = sum(w(:,:,1)*spread(wo(:,i,2)*w(:,i,2),2,ng),1) ! p_12^2
        where (m1.GT.tol) 
          m12(k:l)  = sum(w(:,:,1)*spread(w(:,i,2)*y,2,ng),1) / m1
          co(k:l,1) = m2 / m1
          co(k:l,2) = m0 / m1
        elsewhere
          m12(k:l)    = 0.0
        endwhere
        where (abs(m2).GT.tol) 
          md12(k:l,1) = sum(w(:,:,1)*wo(:,:,1)*
     *                      spread(w(:,i,2)*y,2,ng),1) / m2
          co(k:l,3) = sum(w(:,:,1)*spread(w(:,i,2),2,ng)
     *                        *(wo(:,:,1)**2.),1) / m2
          co(k:l,4) = sum(w(:,:,1)*wo(:,:,1)*
     *                      spread(w(:,i,2)*wo(:,i,2),2,ng),1) / m2
        elsewhere
          md12(k:l,1) = 0.0
        endwhere
        where (abs(m0).GT.tol) 
          md12(k:l,2) = sum(w(:,:,1)*spread(wo(:,i,2)*
     *                        w(:,i,2)*y,2,ng),1) / m0
          co(k:l,5) = sum(w(:,:,1)*spread(w(:,i,2)
     *                        *wo(:,i,2)**2.,2,ng),1) / m0
          co(k:l,6) = sum(w(:,:,1)*wo(:,:,1)*
     *                      spread(w(:,i,2)*wo(:,i,2),2,ng),1) / m0
        elsewhere
          md12(k:l,2) = 0.0
        endwhere
       do jk = 3,d   ! further terms for nuissance directions
        jl = jk-2
        do j=1,ng,1
          where (m1.GT.tol)
            d11(k:l,j,jl) = sum(w(:,:,1)* 
     *                         spread(w(:,i,2)*w(:,j,jk),2,ng),1)
            d13(k:l,j,jl) = sum(w(:,:,1)*wo(:,:,1)*
     *                         spread(w(:,i,2)*w(:,j,jk),2,ng),1)
            d15(k:l,j,jl) = sum(w(:,:,1)*spread(w(:,i,2)*w(:,j,jk)
     *                            *wo(:,i,2),2,ng),1)
            d12(k:l,j,jl) = sum(w(:,:,1)*spread(wo(:,j,jk)*
     *                         w(:,i,2)*w(:,j,jk),2,ng),1)
            d14(k:l,j,jl) = sum(w(:,:,1)*wo(:,:,1)*spread(w(:,i,2)
     *                            *w(:,j,jk)*wo(:,j,jk),2,ng),1)
            d16(k:l,j,jl) = sum(w(:,:,1)*spread(wo(:,i,2)*
     *                         w(:,i,2)*w(:,j,jk)*wo(:,j,jk),2,ng),1)
          end where
        enddo
        where (d23(:,i,:,jl).GT.tol)
          d21(:,i,:,jl) = d11(k:l,:,jl)/d23(:,i,:,jl)
          d22(:,i,:,jl) = d13(k:l,:,jl)/d23(:,i,:,jl)
          d23(:,i,:,jl) = d15(k:l,:,jl)/d23(:,i,:,jl)
        endwhere
        where (abs(d26(:,i,:,jl)).GT.tol)
          d24(:,i,:,jl) = d12(k:l,:,jl)/d26(:,i,:,jl)
          d25(:,i,:,jl) = d14(k:l,:,jl)/d26(:,i,:,jl)
          d26(:,i,:,jl) = d16(k:l,:,jl)/d26(:,i,:,jl)
        endwhere
        where ( (spread(m1,2,ng)) .GT.tol)
          d11(k:l,:,jl) = d11(k:l,:,jl) / spread(m1,2,ng)
          d12(k:l,:,jl) = d12(k:l,:,jl) / spread(m1,2,ng)
        endwhere
        where ( abs(spread(m2,2,ng)) .GT.tol)
          d13(k:l,:,jl) = d13(k:l,:,jl) / spread(m2,2,ng)
          d14(k:l,:,jl) = d14(k:l,:,jl) / spread(m2,2,ng)
        endwhere
        where ( abs(spread(m0,2,ng)) .GT.tol)
          d15(k:l,:,jl) = d15(k:l,:,jl) / spread(m0,2,ng)
          d16(k:l,:,jl) = d16(k:l,:,jl) / spread(m0,2,ng)
        endwhere
       end do ! jl
      enddo
      DEALLOCATE (w)
      df(:,2,:) = df(:,1,:)-df(:,2,:)
      co(:,3) = co(:,1)-co(:,3)
      co(:,6) = co(:,1)-co(:,6)
      where (abs(co(:,3)).GT.tol)
        co(:,4) = (co(:,2)-co(:,4)) / co(:,3)
      elsewhere
        co(:,4) = 0.0
      endwhere
      where (abs(co(:,6)).GT.tol)
        co(:,5) = (co(:,2)-co(:,5)) / co(:,6)
      elsewhere
        co(:,5) = 0.0
      endwhere
      ALLOCATE (e12(gn),e3(ng,dd),ed3(ng,dd),mde(gn,2))
C      Now estimation of m_3 to m_d and m_12 with smooth Backfitting
      e12 = 0.0
      e3  = m3
      ed3 = 0.0
      mde = 0.0
      crit = 1.0
      it   = 0
      m0   = 0.0				! criterion looks at change in m3
      do while ((crit.GT.concr).AND.(maxit.GT.it)) 
       it = it + 1        
       e12  = m12    ! for 2-dim interaction fctn	m12
       mde  = md12   ! for its marginal derivatives
       do j=3,d  
        jl = j-2  
        e12 = e12 - s * (  
     *             matmul( d11(:,2:(ng-1),jl),e3(2:(ng-1),jl) ) +
     *        .5*(e3(1,jl)*d11(:,1,jl)+e3(ng,jl)*d11(:,ng,jl))
     *           + matmul(d12(:,2:(ng-1),jl),ed3(2:(ng-1),jl)) + 
     *        .5*(ed3(1,jl)*d12(:,1,jl)+ed3(ng,jl)*d12(:,ng,jl)) ) 
        mde(:,1) = mde(:,1) - s * (   
     *             matmul(d13(:,2:(ng-1),jl),e3(2:(ng-1),jl)) + 
     *        .5*(e3(1,jl)*d13(:,1,jl)+e3(ng,jl)*d13(:,ng,jl))
     *             + matmul(d14(:,2:(ng-1),jl),ed3(2:(ng-1),jl)) + 
     *        .5*(ed3(1,jl)*d14(:,1,jl)+ed3(ng,jl)*d14(:,ng,jl))  )
        mde(:,2) = mde(:,2) - s * (   
     *             matmul(d15(:,2:(ng-1),jl),e3(2:(ng-1),jl)) + 
     *        .5*(e3(1,jl)*d15(:,1,jl)+e3(ng,jl)*d15(:,ng,jl))
     *             + matmul(d16(:,2:(ng-1),jl),ed3(2:(ng-1),jl)) + 
     *        .5*(ed3(1,jl)*d16(:,1,jl)+ed3(ng,jl)*d16(:,ng,jl))  )
       enddo ! j
 !!!  Auxiliary print outs of the interaction estimates
 !       open(unit=3,file='check.dat',position='APPEND')
 !       write(3,777) (e12(i),i=1,4) 
 !       close(3)	  
        where ( (abs(co(:,3)).GT.tol).AND.
     .          (abs(co(:,6)).GT.tol).AND.
     .          (abs(co(:,4)-co(:,5)).GT.tol) )
          mde(:,2) = (e12-mde(:,1))/co(:,3)-(e12-mde(:,2))/co(:,6)
          mde(:,2) = mde(:,2) / (co(:,4)-co(:,5))
          mde(:,1) = (e12-mde(:,1))/co(:,3)	- mde(:,2)*co(:,4)
          e12  = e12 - mde(:,1)*co(:,1) - mde(:,2)*co(:,2)
        elsewhere
          mde(:,1) = 0.0
          mde(:,2) = 0.0
          e12      = 0.0
        endwhere       
       do j=3,d	! the other additive components (e3) and derivatives (ed3) 
        jk = j-2
        do i=1,ng,1
         	k = 1+(i-1)*ng
          l = i*ng 
          wo(1:ng,:,1) = d21(:,i,:,jk)*spread(e12(k:l),2,ng) +
     *                d22(:,i,:,jk)*spread(mde(k:l,1),2,ng) +
     *                d23(:,i,:,jk)*spread(mde(k:l,2),2,ng)
          wo(i,:,2) = ( sum(wo(2:(ng-1),:,1),1)	+
     *                .5*( wo(1,:,1) + wo(ng,:,1) ) ) *s
          wo(1:ng,:,1) = d24(:,i,:,jk)*spread(e12(k:l),2,ng) +
     *                d25(:,i,:,jk)*spread(mde(k:l,1),2,ng) +
     *                d26(:,i,:,jk)*spread(mde(k:l,2),2,ng)
          wo(i,:,3) = ( sum(wo(2:(ng-1),:,1),1)	+
     *                .5*( wo(1,:,1) + wo(ng,:,1) ) ) *s
        enddo
        e3(:,jk)  = m3(:,jk) - s * (  
     *        sum(wo(2:(ng-1),:,2),1) + .5*(wo(1,:,2)+wo(ng,:,2)) )
        ed3(:,jk) = md3(:,jk) - s * (    
     *        sum(wo(2:(ng-1),:,3),1) + .5*(wo(1,:,3)+wo(ng,:,3)) )
        l = 1
        do k=3,d
         jl = k-2
         if (j.NE.k) then
           e3(:,jk) = e3(:,jk) - s * (
     *             matmul( d31(:,2:(ng-1),jk,l),e3(2:(ng-1),jl) ) +
     *        .5*(e3(1,jl)*d31(:,1,jk,l)+e3(ng,jl)*d31(:,ng,jk,l))
     *           + matmul(d32(:,2:(ng-1),jk,l),ed3(2:(ng-1),jl)) + 
     *        .5*(ed3(1,jl)*d32(:,1,jk,l)+ed3(ng,jl)*d32(:,ng,jk,l)) )
           ed3(:,jk) = ed3(:,jk) - s * (   
     *             matmul(d33(:,2:(ng-1),jk,l),e3(2:(ng-1),jl)) + 
     *        .5*(e3(1,jl)*d33(:,1,jk,l)+e3(ng,jl)*d33(:,ng,jk,l))
     *             + matmul(d34(:,2:(ng-1),jk,l),ed3(2:(ng-1),jl)) + 
     *        .5*(ed3(1,jl)*d34(:,1,jk,l)+ed3(ng,jl)*d34(:,ng,jk,l))  )
           l = l+1
         endif
        enddo ! i,l,jk
        where (abs(df(:,2,jk)).GT.tol)
          ed3(:,jk) = ( e3(:,jk) - ed3(:,jk) ) / df(:,2,jk)
          e3(:,jk)  =  e3(:,jk) - ed3(:,jk) * df(:,1,jk)
        elsewhere
          e3(:,jk)  = 0.0
          ed3(:,jk) = 0.0
        endwhere
       enddo ! j
       crit = sum((sum(e3,2)-m0)**2.,1) / (sum(m0**2.,1)+0.0001) 
       m0   = sum(e3,2) 
      end do
      DEALLOCATE (mde,ed3)
C     Decomposition into m_1, m_2 : start with d12, later test
      ALLOCATE (w(n,ng,2),mde(gn,3))
      do j=1,2,1	           
        wo(:,:,j) = spread(x(:,j),2,ng)-spread(xg(:),1,n)
        w(:,:,j) = (wo(:,:,j) / hd(j)) ** 2.
        where (w(:,:,j).LT.1.0)  
          w(:,:,j) = 0.9375*((1.-w(:,:,j))**2.) / hd(j)
        elsewhere
          w(:,:,j) = 0.0
        endwhere
        y = s*(sum(w(:,2:(ng-1),j),2)+0.5*(w(:,1,j)+w(:,ng,j)))
        do k=1,ng	 ! now projection kernel (must integrate to one)
         where (y.GT.tol)  
           w(:,k,j) = w(:,k,j) / y 
         elsewhere
           w(:,k,j) = 0.0
         endwhere
        enddo
      enddo
      mde = 0.0
      do i=1,ng,1
        k = 1+(i-1)*ng
        l = i*ng
        m1 = sum(w(:,:,1)*spread(w(:,i,2),2,ng),1)  ! p_12
        where (abs(xg).LT.trim)    ! trimming in test statistic 
            mde(k:l,3) = m1 / dn   ! p_12 only used for test statistic
        endwhere        
        where (m1.GT.tol) 
          mde(k:l,1) = m1 / sum(w(:,:,1),1)
          mde(k:l,2) = m1 / sum(w(:,i,2),1)
        endwhere
      enddo
      DEALLOCATE (w)
      crit = 1.0
      it   = 0
      m2   = 0.0
      m0   = 0.0				  ! criterion is change in m1
      do while ((crit.GT.concr).AND.(maxit.GT.it)) 
        it = it + 1
        m1 = 0.0
        do i=2,(ng-1),1
          k = 1+ (i-1)*ng
          l = i*ng
          m1  = m1 + mde(k:l,1)*(e12(k:l)-m2(i))
        enddo
        m1 = m1 + .5*(mde(1:ng,1)*(e12(1:ng)-m2(1)))
     *        + .5*(mde((gn-ng+1):gn,1)*(e12((gn-ng+1):gn)-m2(ng)))
        m1 = m1 *s
        do i=1,ng,1
          k = 2+ (i-1)*ng
          l = i*ng -1
          m2(i) = sum( (e12(k:l)-m1(2:ng-1))*mde(k:l,2) ,1) +
     *          .5*( (e12(k-1)-m1(1)) * mde(k-1,2)
     *               + (e12(l+1)-m1(ng)) * mde(l+1,2) )  
        enddo
        m2 = m2 * s
        crit = sum((m1-m0)**2.,1) / (sum(m0**2.,1)+0.0001) 
        m0   = m1
      end do
C	  Test statistic for Original sample
      do i=1,ng
        k = (i-1)*ng + 2
        j = i*ng -1         
        m0(i) = ( sum( mde(k:j,3) * ( e12(k:j)-
     *          spread(m2(i),1,ng-2) - m1(2:(ng-1)) )**2. ,1)
     *       + .5*( mde(k-1,3)*(e12(k-1)-m2(i)-m1(1) )**2.
     *       +   mde(j+1,3)*(e12(j+1)-m2(i)-m1(ng))**2.) )*s
      enddo
      to =( sum( m0(2:(ng-1)) ,1)+.5*(m0(1)+m0(ng)) ) *s          
      DEALLOCATE (mde,e12)
      ALLOCATE (mde(n,d))
C       Interpolation on real data for m_1, m_2, m_3, ..., m_d
      do j=1,d,1	           
        wo(:,:,j) = spread(x(:,j),2,ng)-spread(xg(:),1,n)
        wo(:,:,j) = wo(:,:,j) / s 
      enddo
      y = 0.0
      mde(1:ng,3:d) = e3
      mde(1:ng,2) = m2
      mde(1:ng,1) = m1      
      do j=1,d
        m0(1:(ng-1)) = mde(2:ng,j) - mde(1:(ng-1),j)        
        m0(ng)       = 0.0
        where ((wo(:,:,j).GE.0.0).AND.(wo(:,:,j).LT.1.0))
	    wo(:,:,j) = wo(:,:,j)*spread(m0,1,n) + spread(mde(1:ng,j),1,n)
        elsewhere
	 wo(:,:,j) = 0.0
        endwhere
        wo(:,1,j) = sum(wo(:,:,j),2)
        where (x(:,j).LT.xg(1))   wo(:,1,j) = mde(1,j)
        where (x(:,j).GT.xg(ng))  wo(:,1,j) = mde(ng,j)
        mde(:,j) = wo(:,1,j)
      enddo
  !!! Auxiliary printout
    !    open(unit=3,file='check.dat')   ! ,position='APPEND')
    !    write(3,777) ((mde(i,j),j=1,d),i=1,n)
    !    close(3)	
      y = sum(mde(:,1:d),2)	! predictor	for H_0 mdoel without constant 
      DEALLOCATE(mde,e3)
      return
      end
CCCCCCCCCCCCCCCCCCCC END CCCCCCCCCCCCCCCCCCCCCCCCCC 