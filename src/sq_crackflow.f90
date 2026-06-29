      subroutine sq_crackflow
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this surboutine modifies surface runoff to account for crack flow

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hhqday(:)   |mm H2O        |surface runoff for the hour in HRUS
 
!!    surfq(:)    |mm H2O        |surface runoff in the HRU for the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hhqday(:)   |mm H2O        |surface runoff for the hour in HRU
!!    surfq(:)    |mm H2O        |surface runoff in the HRU for the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use basin_module
      use hru_module, only : surfq, hhqday, ihru, voltot 
      use soil_module
      use time_module
      
      implicit none

      integer :: j = 0  !none          |HRU number
      real :: voli = 0. !none          |volume available for crack flow
      real :: volcr_eff = 0. !mm       |hydraulically effective crack volume under frozen soil
      real :: frz_prof = 0. !none      |profile hydraulic frozen state
      integer :: ii = 0 !none          |counter

      j = ihru
      if (bsn_cc%froz_soil == 0) then
        volcr_eff = voltot
      else
        frz_prof = Max(0.0, Min(1.0, soil(j)%frz_state)) ** bsn_prm%frz_prof_exp
        volcr_eff = voltot * Max(0.0, Min(1.0, 1.0 - frz_prof))
      end if

      !! subtract hydraulically effective crack flow from surface runoff
      if (surfq(j) > volcr_eff) then
        surfq(j) = surfq(j) - volcr_eff
      else
        surfq(j) = 0.
      endif
      !if (j == 1662) then
      !    print *, "sq_crackflow, voltot: ", voltot, ", surfq:", surfq(j)
      !endif

      if (time%step > 1) then
        voli = 0.
        voli = volcr_eff
        do ii = 1, time%step  !j.jeong 4/24/2009
          if (hhqday(j,ii) > voli) then
            hhqday(j,ii) = hhqday(j,ii) - voli
            voli = 0.
          else
            voli = voli - hhqday(j,ii)
            hhqday(j,ii) = 0.
          endif
        end do
      end if

      return
      end subroutine sq_crackflow