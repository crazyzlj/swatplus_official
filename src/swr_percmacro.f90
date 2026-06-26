      subroutine swr_percmacro
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this surboutine computes percolation by crack flow

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    volcrmin    |mm            |minimum soil volume in profile
!!    voltot      |mm            |total volume of cracks expressed as depth
!!                               |per unit area
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sepbtm(:)   |mm H2O        |percolation from bottom of soil profile for
!!                               |the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Min

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use hru_module, only : sepbtm, voltot, inflpcp, ihru, sepcrk, sepcrktot, volcrmin
      use soil_module
      use basin_module, only : bsn_cc, bsn_prm
      
      implicit none

      external :: layersplit
      integer :: j = 0          !none          |HRU number
      integer :: ly = 0         !none          |counter (soil layer)
      real :: crklch = 0.5      !none          | 
      real :: xx = 0.           !mm H2O        |water deficiency in soil layer
      real :: crk = 0.          !mm H2O        |percolation due to crack flow
      real :: frz_hyd = 0.      !none          |hydraulic activity remaining under frozen soil
      real :: frz_prof = 0.     !none          |profile hydraulic frozen state
      real :: volcr_eff = 0.    !mm            |hydraulically effective crack volume

      j = ihru

      if (bsn_cc%froz_soil == 0) then
          frz_hyd = 1.0
      else
          frz_prof = Max(0.0, Min(1.0, soil(j)%frz_state)) ** bsn_prm%frz_prof_exp
          frz_hyd = Max(0.0, Min(1.0, 1.0 - frz_prof))
      end if
      volcr_eff = voltot * frz_hyd
      sepcrk = Min(volcr_eff, inflpcp)
      sepcrktot = sepcrk
      if (sepcrk > 1.e-4) then
        do ly = soil(j)%nly, 1, -1
          crk = 0.
          xx = 0.
          if (ly == soil(j)%nly) then
          crk = crklch*(soil(j)%ly(ly)%volcr/(soil(j)%phys(ly)%d -          &
                soil(j)%phys(ly-1)%d) * volcr_eff - volcrmin)
          crk = Max(0., crk)
            if (crk < sepcrk) then
              sepcrk = sepcrk - crk
              sepbtm(j) = sepbtm(j) + crk
              soil(j)%ly(ly)%prk = soil(j)%ly(ly)%prk + crk
            else
              sepbtm(j) = sepbtm(j) + sepcrk
              soil(j)%ly(ly)%prk = soil(j)%ly(ly)%prk + sepcrk
              sepcrk = 0.
            end if
          endif
          xx = soil(j)%phys(ly)%fc - soil(j)%phys(ly)%st
          if (xx > 0.) then
            crk = Min(sepcrk, xx)
            soil(j)%phys(ly)%st = soil(j)%phys(ly)%st + crk
            sepcrk = sepcrk - crk
            if (ly /= 1) soil(j)%ly(ly-1)%prk = soil(j)%ly(ly-1)%prk + crk
          end if
          if (sepcrk < 1.e-6) exit
        end do

        !! if soil layers filled and there is still water attributed to
        !! crack flow, it is assumed to percolate out of bottom of profile
        if (sepcrk > 1.e-4) then
          sepbtm(j) = sepbtm(j) + sepcrk
          soil(j)%ly(soil(j)%nly)%prk =                              &                             
                     soil(j)%ly(soil(j)%nly)%prk + sepcrk
        end if
      end if

      return
      end subroutine swr_percmacro