MODULE IPE_Constants_Dictionary

  USE IPE_Precision

  IMPLICIT NONE


  REAL(prec), PARAMETER :: pi      = 3.14159265358979323844_prec
  REAL(prec), PARAMETER :: two_pi  = 2.0_prec * pi
  REAL(prec), PARAMETER :: half_pi = 0.5_prec * pi
  REAL(prec), PARAMETER :: rtd     = 180.0_prec / pi   !radian-->degree
  REAL(prec), PARAMETER :: dtr     = pi / 180.0_prec   !degree-->radian

  REAL(prec), PARAMETER :: kBoltz  = 1.38e-23_prec     ! Boltzmann's constant [J K-1]
  REAL(prec), PARAMETER :: GSCON   = 8.314e+03_prec    ! universal gas constant, as in tucan/ctipe
  REAL(prec), PARAMETER :: AMU     = 1.661e-27_prec    ! Atomic Mass Unit [kg]
  REAL(prec), PARAMETER :: H_mass  = 1.0_prec
  REAL(prec), PARAMETER :: He_mass = 4.0_prec
  REAL(prec), PARAMETER :: N_mass  = 14.0_prec
  REAL(prec), PARAMETER :: N2_mass = 28.0_prec
  REAL(prec), PARAMETER :: O_mass  = 16.0_prec
  REAL(prec), PARAMETER :: O2_mass = 32.0_prec
  REAL(prec), PARAMETER :: NO_mass = N_mass + O_mass

  REAL(prec), DIMENSION(6), PARAMETER :: mass_kg  = (/ O_mass, H_mass, N2_mass, O2_mass, He_mass, N_mass /)
  REAL(prec), DIMENSION(3), PARAMETER :: massn_kg = (/ O_mass, O2_mass, N2_mass /) ! Mass for major neutral species: 1:O; 2:O2; 3:N2

  REAL(prec), PARAMETER :: E0            = 0.035_prec
  REAL(prec), PARAMETER :: WIDTH_Maxwell = 0.050_prec

  REAL(prec), PARAMETER :: earth_radius = 6.3712e+06_prec    !.. Earth radius [meter]
  REAL(prec), PARAMETER :: G0           = 9.80665_prec       !.. strength of the Earth's gravity, nominal average value at the Earth's surface (standard gravity) [m s-2]


  REAL(prec), PARAMETER :: mesh_height_max = 782.0e+03_prec

  ! Parameters for controlling the fixed height grid interpolation
  INTEGER, PARAMETER    :: nheights_geo  = 183
  INTEGER, PARAMETER    :: nlon_geo      = 90
  INTEGER, PARAMETER    :: nlat_geo      = 91
  REAL(prec), PARAMETER :: fillValue_geo = -999999.9999_prec

  ! Parameters governing the size of the emaps, cmaps, and djspectra arrays
  INTEGER, PARAMETER :: n_flux_ipe       = 15
  INTEGER, PARAMETER :: n_bands_ipe      = 21
  INTEGER, PARAMETER :: maps_ipe_size(3) = (/ 21, 20, 7 /)

  ! Conversion factors
  REAL(prec), PARAMETER :: km_to_m     = 1.e+03_prec
  REAL(prec), PARAMETER :: m_to_km     = 1.e-03_prec
  REAL(prec), PARAMETER :: cm_3_to_m_3 = 1.e+06_prec

END MODULE IPE_Constants_Dictionary
