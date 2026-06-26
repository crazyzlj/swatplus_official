      subroutine et_act
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates potential plant transpiration for Priestley-
!!    Taylor and Hargreaves ET methods, and potential and actual soil
!!    evaporation. NO3 movement into surface soil layer due to evaporation
!!    is also calculated.


!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    canstor(:)   |mm H2O        |amount of water held in canopy storage
!!    ep_max       |mm H2O        |maximum amount of transpiration (plant et)
!!                                |that can occur on current day in HRU 
!!    esco(:)      |none          |soil evaporation compensation factor
!!    ihru         |none          |HRU number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    canev        |mm H2O        |amount of water evaporated from canopy
!!                                |storage
!!    ep_max       |mm H2O        |maximum amount of transpiration (plant et)
!!                                |that can occur on current day in HRU
!!    es_day       |mm H2O        |actual amount of evaporation (soil et) that
!!                                |occurs on day in HRU
!!    snoev        |mm H2O        |amount of water in snow lost through
!!                                |sublimation on current day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Min, Max

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
 
      use basin_module
      use organic_mineral_mass_module
      use hru_module, only : hru, canstor, ihru, canev, ep_max,  &
         es_day, pet_day, snoev
      use soil_module
      use plant_module
      use climate_module
      use hydrograph_module
      use water_body_module
      use reservoir_data_module
      
      implicit none

      integer :: j = 0           !none          |HRU number
!!    real, parameter :: esd = 500., etco = 0.80, effnup = 0.1
      real :: esd = 0.           !mm            |maximum soil depth from which evaporation
                                 !              |is allowed to occur
      real :: etco = 0.          !              |
      real :: effnup = 0.        !              ! 
      real :: no3up = 0.         !kg N/ha       |amount of nitrate moving upward in profile 
      real :: es_max = 0.        !mm H2O        |maximum amount of evaporation (soil et)
                                 !              |that can occur on current day in HRU
      real :: eos1 = 0.          !none          |variable to hold intermediate calculation
                                 !              |result
      real :: xx = 0.            !none          |variable to hold intermediate calculation 
                                 !              |result
      real :: cej = 0.           !              |
      real :: eaj = 0.           !none          |weighting factor to adjust PET for impact of
                                 !              |plant cover    
      real :: pet = 0.           !mm H2O        |amount of PET remaining after water stored
                                 !              |in canopy is evaporated
      real :: esleft = 0.        !mm H2O        |potential soil evap that is still available
      real :: evzp = 0.          !              |
      real :: eosl = 0.          !mm H2O        |maximum amount of evaporation that can occur
                                 !              |from soil profile
      real :: dep = 0.           !mm            |soil depth from which evaporation will occur
                                 !              |in current soil layer
      real :: evz = 0.           !              | 
      real :: sev = 0.           !mm H2O        |amount of evaporation from soil layer
      real :: sev_st = 0.        !mm H2O        |evaporation / soil water for no3 flux from layer 1 -> 2
      real :: snowpack_swe = 0.  !mm H2O        |solid snow plus retained liquid water
      real :: snow_evap = 0.     !mm H2O        |evaporation/sublimation removed from snowpack
      real :: frz_surf = 0.      !none          |surface hydraulic frozen state
      real :: frz_evap = 1.      !none          |soil evaporation activity under frozen soil
      real :: cover = 0.         !kg/ha         |soil cover
      real :: wetvol_mm = 0.     !mm            |wetland water volume - average depth over hru
      integer :: ly = 0          !none          |counter     
      integer:: ires = 0  !Jaehak 2022
      integer:: ihyd = 0  !Jaehak 2022

      j = ihru
      pet = pet_day
!!    added statements for test of real statement above
      esd = 500.  !soil(j)%zmx
      etco = 0.80
      effnup = 0.05
      ires= hru(j)%dbs%surf_stor !Jaehak 2022


!! evaporate canopy storage first
!! canopy storage is calculated by the model only if the Green & Ampt
!! method is used to calculate surface runoff. The curve number methods
!! take canopy effects into account in the equations. For either of the
!! CN methods, canstor will always equal zero.
      canev = 0.
      pet = pet - canstor(j)
      if (pet < 0.) then
        canstor(j) = -pet
        canev = pet_day
        pet = 0.
        ep_max = 0.
        es_max = 0.
      else
        canev = canstor(j)
        canstor(j) = 0.
      endif

      if (pet > 1.0e-6) then

        !! compute potential plant evap for methods other that Penman-Monteith
        !if (bsn_cc%pet /= 1) then
          if (pcom(j)%lai_sum <= 3.0) then
            ep_max = pcom(j)%lai_sum * pet / 3.
          else
            ep_max = pet
          end if
          if (ep_max < 0.) ep_max = 0.
        !end if

        !! compute potential soil evaporation
        cej = -5.e-5
        eaj = 0.
        es_max = 0.
        eos1 = 0.
        cover = pl_mass(j)%ab_gr_com%m + pl_mass(j)%rsd_tot%m
        snowpack_swe = hru(j)%sno_mm + hru(j)%sno_liq
        if (snowpack_swe >= 0.5) then
          eaj = 0.5
        else
          eaj = Exp(cej * (cover + 0.1))
        end if
        es_max = pet * eaj
        eos1 = pet / (es_max + ep_max + 1.e-10)
        eos1 = es_max * eos1
        es_max = Min(es_max, eos1)
        es_max = Max(es_max, 0.)
        !if (j == 1662) then
        !  write(9003,*)  "et_act, INPUT, : ires", ires, ", canev: ", canev, ", pet: ", pet, ", ep_max:", ep_max,&
        !      ", eaj: ", eaj, ", eos1:", eos1,", es_max: ", es_max, ", wet_flo: ", wet(j)%flo, ", pet_day:", pet_day
        !end if
        if (wet(j)%flo > 0.) then !wetlands water evaporation reduced by canopy Jaehak 2022
        
          if (pcom(j)%lai_sum <= 4.0) then 
            ihyd = wet_dat(ires)%hyd
            es_max = wet_hyd(ihyd)%evrsv * (1.-pcom(j)%lai_sum / 4.) * pet !adapted from Sakaguchi et al. 2014
          else
            es_max = 0.
          endif

        else  
        
          !! make sure maximum plant and soil ET doesn't exceed potential ET
          if (pet_day < es_max + ep_max) then
            es_max = pet_day - ep_max
            if (pet < es_max + ep_max) then
              es_max = pet * es_max / (es_max + ep_max)
              ep_max = pet * ep_max / (es_max + ep_max)
            end if
          end if
        end if
        
        !! adjust es_max and ep_max for impervous urban cover
        !es_max = 0.5 * es_max
        !ep_max = 0.5 * ep_max
          
        !! initialize soil evaporation variables
        esleft = es_max
        !if (j == 1662) then
        !  write(9003,*)  "init soil et var, esleft:", esleft
        !end if
        !! compute snow evaporation/sublimation from the snowpack before soil evaporation.
        !! With the revised snow routine, snowpack water includes solid snow and retained liquid water.
        snowpack_swe = hru(j)%sno_mm + hru(j)%sno_liq
        if (snowpack_swe > 1.e-6 .and. esleft > 1.e-9) then
          if (w%tave > 0.) then
            snow_evap = Min(esleft, hru(j)%sno_liq)
            hru(j)%sno_liq = hru(j)%sno_liq - snow_evap
            esleft = esleft - snow_evap
            snoev = snoev + snow_evap
            snow_evap = Min(esleft, hru(j)%sno_mm)
            hru(j)%sno_mm = hru(j)%sno_mm - snow_evap
            esleft = esleft - snow_evap
            snoev = snoev + snow_evap
          else
            snow_evap = Min(esleft, hru(j)%sno_mm)
            hru(j)%sno_mm = hru(j)%sno_mm - snow_evap
            esleft = esleft - snow_evap
            snoev = snoev + snow_evap
            snow_evap = Min(esleft, hru(j)%sno_liq)
            hru(j)%sno_liq = hru(j)%sno_liq - snow_evap
            esleft = esleft - snow_evap
            snoev = snoev + snow_evap
          end if
          if (hru(j)%sno_mm < 1.e-9) hru(j)%sno_mm = 0.
          if (hru(j)%sno_liq < 1.e-9) hru(j)%sno_liq = 0.
        endif
        !if (j == 1662) then
        !  write(9003,*)  "after sublimation, esleft:", esleft, ", w%tave:", w%tave
        !end if
        !! compute evaporation from ponded water
        wet_wat_d(j)%evap = 0.
        if (wet(j)%flo > 0.) then
          wetvol_mm = wet(j)%flo / (10. *  hru(j)%area_ha)    !mm=m3/(10.*ha)
          !! take all soil evap from wetland storage before taking from soil
          if (wetvol_mm >= esleft) then
            wetvol_mm = wetvol_mm - esleft
            wet_wat_d(j)%evap = esleft * (10. *  hru(j)%area_ha)
            esleft = 0.
          else
            esleft = esleft - wetvol_mm
            wet_wat_d(j)%evap = wetvol_mm * (10. *  hru(j)%area_ha)
            wetvol_mm = 0.
          endif
          wet(j)%flo = 10. * wetvol_mm * hru(j)%area_ha
          hru(j)%water_evap = wet_wat_d(j)%evap / (10. * hru(j)%area_ha)  !mm=m3/(10*ha)
        endif
        !if (j == 1662) then
        !  write(9003,*)  "after et from pond, esleft:", esleft, ", wet(j)%flo:", wet(j)%flo
        !end if

!! take soil evap from each soil layer
      evzp = 0.
      eosl = esleft
      do ly = 1, soil(j)%nly

        !! depth exceeds max depth for soil evap (esd)
        dep = 0.
        if (ly == 1) then
          dep = soil(j)%phys(1)%d
        else
          dep = soil(j)%phys(ly-1)%d
        endif
        
        if (dep < esd) then
          !! calculate evaporation from soil layer
          evz = eosl * soil(j)%phys(ly)%d / (soil(j)%phys(ly)%d +        &
             Exp(2.374 - .00713 * soil(j)%phys(ly)%d))
          !if (j == 1662) then
          !    write(9003,*) "  d:", soil(j)%phys(ly)%d, ", eosl: ", eosl, ", evz: ", evz
          !endif
          sev = evz - evzp * (1. - hru(j)%hyd%esco)
          evzp = evz
          if (soil(j)%phys(ly)%st < soil(j)%phys(ly)%fc) then
            xx =  2.5 * (soil(j)%phys(ly)%st - soil(j)%phys(ly)%fc) /    &
             soil(j)%phys(ly)%fc
            sev = sev * exp(xx)
          end if
          !if (j == 1662) then
          !    write(9003,*) "  xx:", xx, ", sev: ", sev
          !endif
          frz_surf = Max(0.0, Min(1.0, soil(j)%frz_state)) ** 0.7
          frz_evap = Max(0.05, 1.0 - frz_surf)
          sev = sev * frz_evap
          sev = Min(sev, soil(j)%phys(ly)%st * etco)

          if (sev < 0.) sev = 0.
          if (sev > esleft) sev = esleft

          !! adjust soil storage, potential evap
          if (soil(j)%phys(ly)%st > sev) then
            esleft = esleft - sev
            soil(j)%phys(ly)%st = Max(1.e-6, soil(j)%phys(ly)%st - sev)
          else
            esleft = esleft - soil(j)%phys(ly)%st
            sev = soil(j)%phys(ly)%st
            soil(j)%phys(ly)%st = 0.
          endif
        endif

        !! compute no3 flux from layer 2 to 1 by soil evaporation
        if (ly == 2) then
          if (soil(j)%phys(2)%st > 1.e-3) then
            sev_st = sev / (soil(j)%phys(2)%st)
          else
            sev_st = 0.
          end if
          sev_st = amin1 (1., sev_st)
          no3up = effnup * sev_st * soil1(j)%mn(2)%no3
          no3up = Min(no3up, soil1(j)%mn(2)%no3)
          soil1(j)%mn(2)%no3 = max(0.0001,soil1(j)%mn(2)%no3 - no3up)
          soil1(j)%mn(1)%no3 = max(0.0001,soil1(j)%mn(1)%no3 + no3up)
          !if (j == 1662) then
          !  write(9003,*)  "et_act, sev:", sev, ", soil_st: ", soil(j)%phys(2)%st, "sev_st: ", sev_st, &
          !      ", no3up: ", no3up, ", no3: ", soil1(j)%mn(1)%no3, soil1(j)%mn(2)%no3
          !end if
        endif

      end do    !layer loop

      !! update total soil water content
      soil(j)%sw = 0.
      do ly = 1, soil(j)%nly
        soil(j)%sw = soil(j)%sw + soil(j)%phys(ly)%st
      end do

      !! calculate actual amount of evaporation from soil
      es_day = es_max - esleft
      if (es_day < 0.) es_day = 0.

      end if

      return
      end subroutine et_act