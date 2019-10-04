!
      module module_transf
!----------------------------------------------------------------------- 
!
      PRIVATE
      PUBLIC :: transf
      contains
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: transf
! !INTERFACE:
      subroutine transf
!
! !USES:
      use dynamo_module
      use cons_module,only: dt1dts,rcos0s
      IMPLICIT NONE
!     
! !DESCRIPTION: 
!       - include the boundary condition at the equator
!       - transformation from lam_m to lam_0 (regular grid)
!
! !REVISION HISTORY:
! 05.02.04  <Astrid Maute> <include header> 
! 
! EOP
!
! Local:
      integer :: i,ii,k,kk,j,jj,lat,n
      real :: corfac
!
!
! External:
      real,external :: sddot ! in util.F
!
! From sub fieldline_integrals:
!   zigm11 is int[sig_p*d_1^2/D] ds,   i.e. Sigma_(phi phi)/abs(sin Im)
!   zigm22 is int[sig_p*d_2^2/D] ds,   i.e. Sigma_(lam lam)*abs(sin Im)
!   zigmc  is int[sig_p*d_1*d_2/D] ds, i.e. Sigma_c
!   zigm2  is int[sigma_h] ds,         i.e. Sigma_h
!
!   rim1 is int[(sigma_h-sigma_p*d_1*d_2/D)u_e1 + sigma_p*d_1^2/D u_e2] *A(h_r)*
!                B_e3 ds, i.e.  K_(m phi)^D/abs(sin Im)
!   rim2 is int[(sigma_h+sigma_p*d_1*d_2/D)u_e2 - sigma_p*d_2^2/D u_e1] *A(h_r)*
!                B_e3 ds, K_(m lam)^D ( minus in northern hemisphere
!   Change sign of RIM(2) in S. hemisphere to be compatible with transf
! At this point, rim2 is +-K_(m lam)^D
!
! am 10/04 leave this out and check what difference it makes
! but make sure that field-line integral at equator is definite
!
! Equatorial values:
! Assume that quantities primarily dependent on Pedersen conductivity
!   have field-line integrals 1/4 as large as the averages for next-higher
!   field lines; quantities primarily dependent on Hall conductivity
!   have field-line integrals 0.12 as large as the averages for next-higher
!   field lines.  Exact values chosen should not be important for potential
!   calculation, as long as they are physically reasonable and not too
!   different from adjacent values.
! add factor 1/2 for the equatorial values since both hemispheres are
! added together - two physical points but equator is only one point
!
!      j = kmlat/2+1
!      do i = 1,kmlon
!        zigm11(i,j)= .125*(zigm11(i,j-1)+ zigm11(i,j+1))
!        zigm22(i,j)= .125*(zigm22(i,j-1)+ zigm22(i,j+1))
!        zigmc(i,j) = .125*(zigmc(i,j-1) + zigmc(i,j+1))
!        zigm2(i,j) = .06 *(zigm2(i,j-1) + zigm2(i,j+1))
!        rim(i,j,1) = .06 *(rim(i,j-1,1) + rim(i,j+1,1))
!        rim(i,j,2) = .06 *(rim(i,j-1,2) + rim(i,j+1,2))
!      enddo ! i = 1,kmlon
!
! Include the boundary condition at the equator eq.(5.30) in
! Richmond (1995) Ionospheric Electrodynamics use. Mag. Apex Coord. 
!   
! Sigma_(phi phi)/|sin Im| = 0.5*Sigma_cowling/|sin Im|
!        = 0.5/ |sin Im|*[Simag_(phi phi) -
!                        Sigma_(phi lam)*Sigma_(lam phi)/Sigma_(lam lam)]
!        = 0.5/|sin Im|*[Simag_(phi phi) + 
!                       (Sigma_h-Sigma_c)*(Sigma_h + Sigma_c)/Sigma_(lam lam)]
!  K_(m phi)^D / |sin I_m| =  1/2/|sin Im|*(K_(m phi) - 
!                       Sigma_(phi lam)/Sigma_(lam lam)*K_(m lam)^D)
! factor 1/2 is taken care of above when the conductances at the
! equator are set
!
      j = kmlath      ! kmlath = (kmlat+1)/2
      do i = 1,kmlon
        zigm11(i,j) = zigm11(i,j)+ (zigm2(i,j)-zigmc(i,j))*
     |                (zigm2(i,j)+zigmc(i,j))/zigm22(i,j)
        rim(i,j,1)  = rim(i,j,1) - (zigm2(i,j)-zigmc(i,j))/
     |                zigm22(i,j)*rim(i,j,2)
        zigm11(i,j) = 0.5*zigm11(i,j)
        rim(i,j,1)  = 0.5*rim(i,j,1)
      enddo

!        
! Using notation of Richmond (1995) on right-hand side below:
! Sigma_(phi phi)/ |sin I_m|  = zigm11
! Sigma_(lam lam) * |sin I_m| = zigm22
! Sigma_(phi lam)             = +/-(zigm2-zigmc)
! Sigma_(lam phi)             = -/+(zigm2+zigmc)
! K_(m phi)^D / |sin I_m|     = rim(1)
! K_(m lam)^D                 = +/-rim(2)
!
! Transforming PDE from original apex (lam_m) to new apex grid (lam_0)
!     lam_m is irregular spaced in latitude
!     lam_0 is regular spaced in latitude (used for derivatives)
! the whole PDE is divided by d lam_m /d lam_0
! 
! sign of K_(m lam)^D in southern hemisphere remains reversed.
! for the mixed terms the transformation factors cancel out (zigmc, zigm2)
! DT1DTS : d lam_0/ d lam_m / |sin I_m|
! RCOS0S : cos(lam_0)/ cos(lam_m)
!
! corfac: |sin I_m|*d lam_m/d lam_0 * cos(lam_0)/ cos(lam_m)
! zigm11: |sin I_m|*d lam_m/d lam_0 * cos(lam_0)/ cos(lam_m)
! zigm22: 1/|sin I_m|*d lam_0/d lam_m * cos(lam_m)/ cos(lam_0)
! rim(1): |sin I_m|*d lam_m/d lam_0
! rim(2): cos(lam_m)/ cos(lam_0)        
!
      do j=2,kmlat-1
        corfac = rcos0s(j)/dt1dts(j)
        do i=1,kmlon
          zigm11(i,j) = zigm11(i,j)*corfac
          zigm22(i,j) = zigm22(i,j)/corfac
          rim(i,j,1)  = rim(i,j,1)/dt1dts(j) 
          rim(i,j,2)  = rim(i,j,2)/rcos0s(j) 
        enddo ! i,kmlon
      enddo ! j=2,kmlat-1
!
! Periodic points:
      do j=1,kmlat
        zigm11(kmlonp1,j) = zigm11(1,j)
        zigmc (kmlonp1,j) = zigmc (1,j)
        zigm2 (kmlonp1,j) = zigm2 (1,j)
        zigm22(kmlonp1,j) = zigm22(1,j)
        rim(kmlonp1,j,:)  = rim(1,j,:)
      enddo ! j=1,kmlat
!
! zigm11 = Sigma_(phi phi)(0)= Sigma_(phi phi)(m) * cos lam_0 /cos lam_m * d lam_m /d lam_0
! zigm22 = Sigma_(lam lam)(0)= Sigma_(lam lam)(m) * cos lam_m /cos lam_0 * d lam_0 /d lam_m
! +-(zigm2-zigmc)= Sigma_(phi lam)(0) = Sigma_(phi lam)(m)
! -+(zigm2+zigmc)= Sigma_(lam phi)(0) = Sigma_(lam phi)(m)
! rim(1) = K_(m phi)^D(0)    = K_(m phi)^D(m) * d lam_m /d lam_0   
! rim(2) = +/-K_(m lam)^D(0) = +/-K_(m lam)^D(m) * cos lam_m /cos lam_0 
!  
! Compute polar values for the conductances, 4th order interpolation:
! 
      zigm11(1,    1) = (4.*sddot(kmlon,unitvm,zigm11(1,      2))-
     1  sddot(kmlon,unitvm,zigm11(1,      3)))/(3.*float(kmlon))
      zigm11(1,kmlat) = (4.*sddot(kmlon,unitvm,zigm11(1,kmlat-1))-
     1  sddot(kmlon,unitvm,zigm11(1,kmlat-2)))/(3.*float(kmlon))
      zigmc(1,    1) = (4.*sddot(kmlon,unitvm,zigmc(1,      2))-
     1  sddot(kmlon,unitvm,zigmc(1,      3)))/(3.*float(kmlon))
      zigmc(1,kmlat) = (4.*sddot(kmlon,unitvm,zigmc(1,kmlat-1))-
     1  sddot(kmlon,unitvm,zigmc(1,kmlat-2)))/(3.*float(kmlon))
      zigm2(1,    1) = (4.*sddot(kmlon,unitvm,zigm2(1,      2))-
     1  sddot(kmlon,unitvm,zigm2(1,      3)))/(3.*float(kmlon))
      zigm2(1,kmlat) = (4.*sddot(kmlon,unitvm,zigm2(1,kmlat-1))-
     1  sddot(kmlon,unitvm,zigm2(1,kmlat-2)))/(3.*float(kmlon))
      zigm22(1,    1) = (4.*sddot(kmlon,unitvm,zigm22(1,      2))-
     1  sddot(kmlon,unitvm,zigm22(1,      3)))/(3.*float(kmlon))
      zigm22(1,kmlat) = (4.*sddot(kmlon,unitvm,zigm22(1,kmlat-1))-
     1  sddot(kmlon,unitvm,zigm22(1,kmlat-2)))/(3.*float(kmlon))
! 
! Extend over longitude                                        
      do i = 2,kmlon
        zigm11(i,    1)  = zigm11(1,    1)
        zigm11(i,kmlat)  = zigm11(1,kmlat)
        zigmc(i,     1)  = zigmc(1,     1)
        zigmc(i, kmlat)  = zigmc(1, kmlat)
        zigm2(i,     1)  = zigm2(1,     1)
        zigm2(i, kmlat)  = zigm2(1, kmlat)
        zigm22(i,     1) = zigm22(1,    1)
        zigm22(i, kmlat) = zigm22(1,kmlat)
      enddo ! i = 2,kmlon
!
! RHS vector (I_1,I_2): average over poles:
      do i = 1,kmlon
        rim(i,1,1) = .5*(rim(i,2,1)-rim(1+mod(i-1+kmlon/2,kmlon),2,1))
        rim(i,kmlat,1) = .5*(rim(i,kmlat-1,1)-
     |    rim(1+mod(i-1+kmlon/2,kmlon),kmlat-1,1))
        rim(i,1,2) = .5*(rim(i,2,2)-rim(1+mod(i-1+kmlon/2,kmlon),2,2))
        rim(i,kmlat,2) = .5*(rim(i,kmlat-1,2)-
     |    rim(1+mod(i-1+kmlon/2,kmlon),kmlat-1,2))
      enddo ! i = 1,kmlon
! 
! Periodic points:
      do j=1,kmlat
        zigm11(kmlonp1,j) = zigm11(1,j)
        zigmc (kmlonp1,j) = zigmc (1,j)
        zigm2 (kmlonp1,j) = zigm2 (1,j)
        zigm22(kmlonp1,j) = zigm22(1,j)
        rim(kmlonp1,j,:)  = rim(1,j,:)
      enddo ! j=1,kmlat
!
! Save to secondary history:
!     real,dimension(kmlonp1,kmlat) ::
!    |  zigm11,  ! sigma11*cos(theta0)
!    |  zigmc,   ! sigmac
!    |  zigm2,   ! sigma2
!    |  zigm22   ! sigma22/cos(theta0)
!
!     do j=1,kmlat
!       do i=1,kmlonp1
!         zigm11_plt(i,:) = zigm11(i,j)
!         zigmc_plt (i,:) = zigmc (i,j)
!         zigm2_plt (i,:) = zigm2 (i,j)
!         zigm22_plt(i,:) = zigm22(i,j)
!         rim1_plt  (i,:) = rim   (i,j,1)
!         rim2_plt  (i,:) = rim   (i,j,2)
!       enddo ! i=1,kmlonp1
! 
! Folding south on to northern hemisphere and calculation of RHS has 
!   been moved to subroutine dynamo.
!
      end subroutine transf
!-----------------------------------------------------------------------
      end module module_transf
!-----------------------------------------------------------------------
