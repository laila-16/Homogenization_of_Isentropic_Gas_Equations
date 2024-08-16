! =====================================================
subroutine rp1(maxm,meqn,mwaves,maux,mbc,mx,zl,zr,auxl,auxr,fwave,s,amdz,apdz)
! =====================================================

!     # Riemann solver for the isentropic gas equations in 1d, with
!     #  variable coefficients (heterogeneous media)

! waves:     2
! equations: 2
! aux fields: 1

! Conserved quantities:
!       1 density
!       2 flux

! Auxiliary variables:
!       1 cross_section
!

!     # On input, zl contains the state vector at the left edge of each cell
!     #           zr contains the state vector at the right edge of each cell

!     # On output, wave contains the waves,
!     #            s the speeds,
!     #
!     #            amdq = A^- Delta z,
!     #            apdq = A^+ Delta z,
!     #                   the decomposition of the flux difference
!     #                       f(zr(i-1)) - f(zl(i))
!     #                   into leftgoing and rightgoing parts respectively.
!     #

!     # Note that the i'th Riemann problem has left state zr(:,i-1)
!     #                                    and right state zl(:,i)
!     # From the basic clawpack routines, this routine is called with zl = zr

    implicit none 
    
    integer, intent(in) :: maxm, meqn, mwaves, mbc, maux, mx
    real(kind=8), intent(in) :: zl(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: zr(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: auxl(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: auxr(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(out) :: s(mwaves, 1-mbc:maxm+mbc)
    real(kind=8), intent(out) :: fwave(meqn, mwaves, 1-mbc:maxm+mbc)
    real(kind=8), intent(out) :: amdz(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(out) :: apdz(meqn, 1-mbc:maxm+mbc)
    


!     local arrays
!     ------------
    real(kind=8) :: fluxdiff(2), RM(2,2), RM_inv(2,2)
    real(kind=8) :: lambdabar(2) !Eigenvalues to compute speeds
    real(kind=8) :: betaM(2) !Coefficients of eigenvectors to compute waves
    integer :: i, m
    real(kind=8) :: rhol, rhom, rhor, rhoul, rhoum, rhour , al, am, ar,p_primel, p_primem, p_primer, pl, pm, pr


! Parameters
    double precision :: kappa, gamma

    kappa = 1.d0
    gamma = 1.4d0

    do i = 2-mbc, mx+mbc
        !Vector of consereved quantities
        rhor = zl(1,i)
        rhour = zl(2,i)
        rhol = zr(1,i-1)
        rhoul = zr(2,i-1)

        !Spatial varying function
        ar = auxl(1,i)
        al = auxr(1,i-1)

        !p prime
        pL = kappa*rhol**(gamma)
        pR = kappa*rhor**(gamma)
        p_primeL = kappa*(gamma)*rhol**(gamma-1)
        p_primeR = kappa*(gamma)*rhor**(gamma-1)

        !Averages
        rhom = 0.5d0*(rhor+rhol)
        rhoum = 0.5d0*(rhour+rhoul)
        am = 0.5d0*(al + ar)
        pM = 0.5d0*(pL + pR)
        p_primeM = 0.5d0*(p_primeL + p_primeR)


        !Eigenvalues
        lambdabar(1) = min((al*rhoul/rhol)-al*sqrt(p_primeL), (ar*rhour/rhor)-ar*sqrt(p_primeR))
        lambdabar(2) = max((al*rhoul/rhol)+al*sqrt(p_primeL), (ar*rhour/rhor)+ar*sqrt(p_primeR))
        
   
        !flux difference
        fluxdiff(1) = ar*rhour - al*rhoul
        fluxdiff(2) = ar*rhour**2/rhor + ar*pR - al*rhoul**2/rhol - al*pL - pM*(ar -al)
  
        !Right eigenvectors
        RM(1,1) = rhom
        RM(1,2) = rhom
        RM(2,1) = rhoum - rhom*sqrt(p_primeM)
        RM(2,2) = rhoum + rhom*sqrt(p_primeM)
        
        !Inverse of right eigenvectors
        RM_inv(1,1) = 0.5d0*(rhoum + rhom*sqrt(p_primeM))/(sqrt(p_primeM)*rhom**2)
        RM_inv(2,1) = 0.5d0*(-rhoum + rhom*sqrt(p_primeM))/(sqrt(p_primeM)*rhom**2)
        RM_inv(1,2) = 0.5d0*(-1./(rhom*sqrt(p_primeM)))
        RM_inv(2,2) = 0.5d0*(1./(rhom*sqrt(p_primeM)))
        
	!Coefficient of eigenvectors
        betaM(1) = RM_inv(1, 1)*fluxdiff(1)+ RM_inv(1, 2)*fluxdiff(2)
        betaM(2) = RM_inv(2, 1)*fluxdiff(1)+ RM_inv(2, 2)*fluxdiff(2)

        !Flux waves
        fwave(:,1,i) = betaM(1)*RM(:,1)
        fwave(:,2,i) = betaM(2)*RM(:,2)

        !Speeds
        s(1,i) = lambdabar(1)
        s(2,i) = lambdabar(2)
    
    END DO

!     # compute the leftgoing and rightgoing fluctuations:
!     # Note s(1,i) < 0   and   s(2,i) > 0.
    apdz(:,:) = 0.d0
    amdz(:,:) = 0.d0

    do m=1,meqn
        do i = 2-mbc, mx+mbc
            if (s(1,i)>0) then
                apdz(m,i) = apdz(m,i)+ fwave(m,1,i)
            else
                amdz(m,i) = amdz(m,i) + fwave(m,1,i)
            end if 
            if (s(2,i)<0) then
                amdz(m,i) = amdz(m,i) + fwave(m,2,i)
            else
                apdz(m,i) = apdz(m,i)+ fwave(m,2,i)
            end if
        end do
    end do
    return
    end subroutine rp1
