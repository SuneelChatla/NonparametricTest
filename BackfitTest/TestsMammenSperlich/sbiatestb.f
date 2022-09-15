C  Stefan Sperlich, created in 2018  (from mamt2.f subroutine backfb)
C
C  this file contains only the subroutine for the bootstrap part
C   
C     All X are on similar scale around interval [-ra/2;ra/2]
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      SUBROUTINE  sbiatestb   	for loc.lin. smooth backfit
C							    for bootstrap samples 
C                is like backia, loc.lin. smooth backfit
C                with interaction m_12 and its decomposition
C                Calculates test statistic value
C  INPUTS :
C   x   = designmatrix, 
C   y   = dependent observations
C   yb  = cond expect under H0 (potentially without constant)
C   to  = value of original statistic
C   rd  = random seed for replicability of simulations
C   h   = bandwidth vector for one-dim components
C   hd  = bandwidth vector for interaction m_12 (called e_12)
C   n   = number of observations
C   d   = number of dimensions,    
C   ng  = grid size in each direction
C   ra  = range of X on which integral is calculated
C   trim = threshold for trimming (see papers with Dette and Tjostheim)
C  OUTPUT :
C   pv =  p-value of test
C
      subroutine sbiatestb(x,y,yb,n,d,h,hd,pv,ng,nb,rd,to,trim,ra) 
      implicit real(kind=8) (a-h,o-z)
      integer  n,d,it,i,j,k,l,ng,maxit,jk,jl,nb,lb,dd,gn
      parameter(maxit=15,concr=0.00001,tol=0.00001)
      real(kind=8)  x(n,d),y(n),xg(ng),h(d),s,hd(2),ra
     .    ,d11(ng*ng,ng,d-2),d21(ng,ng,ng,d-2),d31(ng,ng,d-2,d-3)
     .    ,d12(ng*ng,ng,d-2),d22(ng,ng,ng,d-2),d32(ng,ng,d-2,d-3)
     .    ,d13(ng*ng,ng,d-2),d23(ng,ng,ng,d-2),d33(ng,ng,d-2,d-3)
     .    ,d14(ng*ng,ng,d-2),d24(ng,ng,ng,d-2),d34(ng,ng,d-2,d-3)
     .    ,d15(ng*ng,ng,d-2),d25(ng,ng,ng,d-2)
     .    ,d16(ng*ng,ng,d-2),d26(ng,ng,ng,d-2) ,d3(ng*ng,3) ! old notation
     .    ,co(ng*ng,6),df(ng,2,d-2),m0(ng),m1(ng),m2(ng),wo(n,ng,d)
     .    ,m12(ng*ng),md12(ng*ng,2),m3(ng,d-2),md3(ng,d-2)    
     .    ,w(n,ng,d),e12(ng*ng),e3(ng,d-2),ed3(ng,d-2),mde(ng*ng,2)
     .    ,yb(n),mb(n,2),pv,w1,w2,pw,rd,to,tb,trim
        gn = ng*ng
        dd = d-2
      dn = float(n)
      cy = sum(y,1) / dn
      y = y - cy
      s = ra / (ng-1.0)    
      do j=1,ng,1	           
        xg(j) = (j - 1.0)*s  - 0.5*ra  
      enddo
C      creation of terms d3 neede for decomposition & test statistic
       do j=1,2,1	           
        wo(:,:,j) = spread(x(:,j),2,ng)-spread(xg(:),1,n)
        w(:,:,j) = (wo(:,:,j) / hd(j)) ** 2.
        where (w(:,:,j).LT.1.0)  
          w(:,:,j) = 0.9375*((1.-w(:,:,j))**2.) / hd(j)
        elsewhere
          w(:,:,j) = 0.0
        endwhere	   ! now projection kernel (must integrate to one)
        mb(:,1) =s*(sum(w(:,2:(ng-1),j),2)+0.5*(w(:,1,j)+w(:,ng,j)))
        do k=1,ng
         where (mb(:,1).GT.tol)  
           w(:,k,j) = w(:,k,j) / mb(:,1) 
         elsewhere
           w(:,k,j) = 0.0
         endwhere
        enddo
       enddo
       d3 = 0.0
       do i=1,ng,1
         k = 1+(i-1)*ng
         l = i*ng
         m1 = sum(w(:,:,1)*spread(w(:,i,2),2,ng),1)  ! p_12
         where (abs(xg).LT.trim)   ! trimming in test statistic 
            d3(k:l,3) = m1 / dn   ! p_12 only used for test statistic
         endwhere        
         where (m1.GT.tol) 
           d3(k:l,1) = m1 / sum(w(:,:,1),1)
           d3(k:l,2) = m1 / sum(w(:,i,2),1)
         endwhere
       enddo
C      now kernel weights for SB with interaction
      do j=1,d,1	           
        wo(:,:,j) = spread(x(:,j),2,ng)-spread(xg(:),1,n)
        w(:,:,j) = (wo(:,:,j) / h(j)) ** 2.
        where (w(:,:,j).LT.1.0)  
          w(:,:,j) = 0.9375*((1.-w(:,:,j))**2.) / h(j)
        elsewhere
          w(:,:,j) = 0.0
        endwhere	   ! now projection kernel (must integrate to one)
        mb(:,1) =s*(sum(w(:,2:(ng-1),j),2)+0.5*(w(:,1,j)+w(:,ng,j)))
        do k=1,ng
         where (mb(:,1).GT.tol)  
           w(:,k,j) = w(:,k,j) / mb(:,1) 
         elsewhere
           w(:,k,j) = 0.0
         endwhere
        enddo
      enddo
C      auxiliary terms/weights like densities, starting with dimensions>2
      do j=3,d
       it = j-2
       m1 = sum(w(:,:,j),1)              ! for p_3
       m2 = sum(w(:,:,j)*wo(:,:,j),1)    ! for p_3^3  hilfsf
       where (m1.GT.tol)
         df(:,1,it) = m2 / m1
       elsewhere
         df(:,1,it) = 0.0	 !
       endwhere
       where (abs(m2).GT.tol)
         df(:,2,it) = sum(w(:,:,j)*wo(:,:,j)**2.,1) / m2
       elsewhere
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
      d1 = 0.0      
      do i=1,ng,1
        k = 1+(i-1)*ng
        l = i*ng
        m1 = sum(w(:,:,1)*spread(w(:,i,2),2,ng),1)  ! p_12
        m2 = sum(w(:,:,1)*wo(:,:,1)*spread(w(:,i,2),2,ng),1) ! p_12^1
        m0 = sum(w(:,:,1)*spread(wo(:,i,2)*w(:,i,2),2,ng),1) ! p_12^2
        where (m1.GT.tol) 
          co(k:l,1) = m2 / m1
          co(k:l,2) = m0 / m1
        endwhere
        where (abs(m2).GT.tol) 
          co(k:l,3) = sum(w(:,:,1)*spread(w(:,i,2),2,ng)
     *                        *(wo(:,:,1)**2.),1) / m2
          co(k:l,4) = sum(w(:,:,1)*wo(:,:,1)*
     *                      spread(w(:,i,2)*wo(:,i,2),2,ng),1) / m2
        endwhere
        where (abs(m0).GT.tol) 
          co(k:l,5) = sum(w(:,:,1)*spread(w(:,i,2)
     *                        *wo(:,i,2)**2.,2,ng),1) / m0
          co(k:l,6) = sum(w(:,:,1)*wo(:,:,1)*
     *                      spread(w(:,i,2)*wo(:,i,2),2,ng),1) / m0
        endwhere
       do jk = 3,d
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
        where ( (spread(m1,2,ng)).GT.tol)
          d11(k:l,:,jl) = d11(k:l,:,jl) / spread(m1,2,ng)
          d12(k:l,:,jl) = d12(k:l,:,jl) / spread(m1,2,ng)
        endwhere
        where ( abs(spread(m2,2,ng)).GT.tol)
          d13(k:l,:,jl) = d13(k:l,:,jl) / spread(m2,2,ng)
          d14(k:l,:,jl) = d14(k:l,:,jl) / spread(m2,2,ng)
        endwhere
        where ( abs(spread(m0,2,ng)).GT.tol)
          d15(k:l,:,jl) = d15(k:l,:,jl) / spread(m0,2,ng)
          d16(k:l,:,jl) = d16(k:l,:,jl) / spread(m0,2,ng)
        endwhere
       end do ! jl
      enddo
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Now the  golden cut Bootstrap
      mb(:,1) = yb	            ! predictor H_0 without constant 
      mb(:,2) = y - mb(:,1)		! residuen under H_0
      pv = 0.0
      w1 = (1.-sqrt(5.))/2.
      w2 = (1.+sqrt(5.))/2.
      pw = (1.+sqrt(5.))/(2.*sqrt(5.))           
      S2P31M = 2147483647.0
      do lb=1,nb		! loop over bootstrap replicates
        SEED = rd	    ! generating uniform random vars
        DO 5 i=1,n
          SEED = MOD(16807.d0*SEED,S2P31M)
 5        yb(i) = SEED / 2.0**31
        rd = SEED 	 ! end of random generator 
        where (yb.le.pw)
          yb = mb(:,1) + mb(:,2)*w1
        elsewhere
          yb = mb(:,1) + mb(:,2)*w2
        endwhere
        yb    = yb - ( sum(yb,1)/n )  ! create centered bootstrap y*
C      calculate e_12, m_1, m_2 for Bootstrap-samples
       do j=1,d,1	           
        wo(:,:,j) = spread(x(:,j),2,ng)-spread(xg,1,n)
       enddo
C      create marginal starting terms (first for d=3) in BOOTSTRAP
       do j=3,d
        it = j-2
        m1 = sum(w(:,:,j),1)			 !  p_3
        m2 = sum(w(:,:,j)*wo(:,:,j),1)   ! p_3^3  hilfsf
        where (m1.GT.tol)
          m3(:,it)  = sum(w(:,:,j)*spread(yb,2,ng),1) / m1
        elsewhere
          m3(:,it)   = 0.0	 !	m_3
        endwhere
        where (abs(m2).GT.tol)
          md3(:,it) = sum(w(:,:,j)*wo(:,:,j)*spread(yb,2,ng),1) / m2
        elsewhere
          md3(:,it)  = 0.0	 !	m_3^3  hilfsf
        endwhere
       enddo ! j
       do i=1,ng,1
         k = 1+(i-1)*ng
         l = i*ng
         m1 = sum(w(:,:,1)*spread(w(:,i,2),2,ng),1)  ! p_12
         m2 = sum(w(:,:,1)*wo(:,:,1)*spread(w(:,i,2),2,ng),1) ! p_12^1
         m0 = sum(w(:,:,1)*spread(wo(:,i,2)*w(:,i,2),2,ng),1) ! p_12^2
         where (m1.GT.tol) 
           m12(k:l)  = sum(w(:,:,1)*spread(w(:,i,2)*yb,2,ng),1) / m1
         elsewhere
           m12(k:l)    = 0.0
         endwhere
         where (abs(m2).GT.tol) 
           md12(k:l,1) = sum(w(:,:,1)*wo(:,:,1)*
     *                      spread(w(:,i,2)*yb,2,ng),1) / m2
         elsewhere
           md12(k:l,1) = 0.0
         endwhere
         where (abs(m0).GT.tol) 
           md12(k:l,2) = sum(w(:,:,1)*spread(wo(:,i,2)*
     *                        w(:,i,2)*yb,2,ng),1) / m0
         elsewhere
           md12(k:l,2) = 0.0
         endwhere
       enddo
C      Now estimate m_3, m_12 with Backfitting	for BOOTSTRAP
       e3   = m3
       ed3  = 0.0
       mde  = 0.0
       crit = 1.0
       it   = 0
       m0   = 0.0				   ! criterion in change in m3
       do while ((crit.GT.concr).AND.(maxit.GT.it)) 
        it = it + 1        
        e12  = m12 
        mde  = md12
        do j=3,d  
         jl = j-2  
         e12 = e12 - s * (  
     *             matmul( d11(:,2:(ng-1),jl),e3(2:(ng-1),jl) ) +
     *        .5*(e3(1,jl)*d11(:,1,jl)+e3(ng,jl)*d11(:,ng,jl))
     *           + matmul(d12(:,2:(ng-1),jl),ed3(2:(ng-1),jl)) + 
     *        .5*(ed3(1,jl)*d12(:,1,jl)+ed3(ng,jl)*d12(:,ng,jl))  )
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
        do j=3,d
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
C     Decomposition of e_12 into m_1, m_2
       crit = 1.0
       it   = 0
       m2   = 0.0
       m0   = 0.0			   ! criterion in change in m1
       do while ((crit.GT.concr).AND.(maxit.GT.it)) 
        it = it + 1
        m1 = 0.0
        do i=2,(ng-1),1
          k = 1+ (i-1)*ng
          l = i*ng
          m1  = m1 + d3(k:l,1)*(e12(k:l)-m2(i))
        enddo
        m1 = m1 + .5*(d3(1:ng,1)*(e12(1:ng)-m2(1)))
     *        + .5*(d3((gn-ng+1):gn,1)*(e12((gn-ng+1):gn)-m2(ng)))
        m1 = m1 *s
        do i=1,ng,1
          k = 2+ (i-1)*ng
          l = i*ng -1
          m2(i) = sum( (e12(k:l)-m1(2:ng-1))*d3(k:l,2),1) +
     *          .5*( (e12(k-1)-m1(1)) * d3(k-1,2)
     *               + (e12(l+1)-m1(ng)) * d3(l+1,2) )  
        enddo
        m2 = m2 * s
        crit = sum((m1-m0)**2.,1) / (sum(m0**2.,1)+0.0001) 
        m0   = m1
       end do
C	  Test statistic for BOOTSTRAP samples
       do i=1,ng
         k = (i-1)*ng + 2
         j = i*ng -1         
         m0(i) = ( sum( d3(k:j,3) * ( e12(k:j)-
     *          spread(m2(i),1,ng-2) - m1(2:(ng-1)) )**2. ,1)
     *       + .5*( d3(k-1,3)*(e12(k-1)-m2(i)-m1(1) )**2.
     *       +   d3(j+1,3)*(e12(j+1)-m2(i)-m1(ng))**2.) )*s
       enddo
       tb =( sum( m0(2:(ng-1)) ,1)+.5*(m0(1)+m0(ng)) ) *s    
       if (tb.GT.to)  pv = pv+1.0  ! compare with original test
      enddo		   ! loop over bootstrap replicates
      pv = pv / float(nb)	 ! calculate bootstrap  p-value
      return
      end
CCCCCCCCCCCCCCCCCCCC END CCCCCCCCCCCCCCCCCCCCCCCCCC 