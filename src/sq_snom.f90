      subroutine sq_snom
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine predicts daily snom melt when the maximum air
!!    temperature exceeds snow melt temperature (0.5 degrees Celsius),
!!    considering rain-on-snow advection and snowpack ripening.

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ihru         |none          |HRU number
!!    snocov1      |none          |1st shape parameter for snow cover equation
!!                                |This parameter is determined by solving the
!!                                |equation for 50% snow cover
!!    snocov2      |none          |2nd shape parameter for snow cover equation
!!                                |This parameter is determined by solving the
!!                                |equation for 95% snow cover
!!    snotmp       |deg C         |temperature of snow pack in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    wst(:)%weat%ts(:)  |mm H2O        |precipitation for the time step during day
!!    snofall      |mm H2O        |amount of precipitation falling as freezing rain/snow on day
!!    snomlt       |mm H2O        |amount of water in snow melt for the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Real, Sin, Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use time_module
      use hydrograph_module
      use hru_module, only : hru, ihru, precip_eff, snofall, snomlt 
      use climate_module, only:  w
      use output_landscape_module
      
      implicit none

      integer :: j = 0      !none       |HRU number
      real :: smfac = 0.    !           |
      real :: rto_sno  = 0. !none       |ratio of current day's snow water to minimum amount needed to
                            !           |cover ground completely 
      real :: snocov = 0.   !none       |fraction of HRU area covered with snow
      real :: snotmp = 0.   !deg C      |temperature of snow pack
      real :: snomlt_ros = 0. !         |
      real :: rain_mm = 0.    !mm       |temporary variable
      real :: c_ros = 0.0125  !none     |Rain-On-Snow advective heat coefficient
      real :: liq_max_ratio = 0.1 !none |max. liquid retention capacity, 10% of the snowpack
      real :: liq_capacity = 0. !       |
      real :: total_liquid_input = 0. ! |

      snotmp = 0.
      snofall = 0.
      snomlt = 0.
      snomlt_ros = 0.
      rain_mm = 0.
      c_ros = 0.0125
      liq_max_ratio = 0.1
      liq_capacity = 0.
      total_liquid_input = 0.

      j = ihru
      
      !! calculate snow fall
      if (w%tave <= hru(j)%sno%falltmp) then
        hru(j)%sno_mm = hru(j)%sno_mm + precip_eff
        snofall = precip_eff
        precip_eff = 0.
        !! set subdaily effective precip to zero
        if (time%step > 1) w%ts = 0.
      endif
      
      !! process existing snowpack
      if (hru(j)%sno_mm > 0.) then
        !! estimate internal snow pack temperature
        snotmp = hru(j)%sno_tmp
        snotmp = snotmp * (1. - hru(j)%sno%timp) + w%tave * hru(j)%sno%timp
        hru(j)%sno_tmp = snotmp
 
        !! calculate snow melt
        if (w%tmax > hru(j)%sno%melttmp) then
          !! adjust degree-day melt factor for time of year
          smfac = (hru(j)%sno%meltmx + hru(j)%sno%meltmn) / 2. + Sin((time%day - 81) / 58.09) *     &
                        (hru(j)%sno%meltmx - hru(j)%sno%meltmn) / 2.        !! 365/2pi = 58.09
          snomlt = smfac * (((snotmp + w%tmax)/2.) - hru(j)%sno%melttmp)
          if (snomlt < 0.) snomlt = 0.
          
          !! calculate rain-on-snow advection heat
          if (w%tave > hru(j)%sno%falltmp) then
            rain_mm = precip_eff
            if (rain_mm > 0. .and. w%tave > 0.) then
              if (snotmp >= -1.0) then
                  snomlt_ros = rain_mm * c_ros * w%tave
              else
                  snomlt_ros = (rain_mm * c_ros * w%tave) + (snotmp * 0.1)
              endif
              if (snomlt_ros < 0.) snomlt_ros = 0.
            endif 
          endif
          
          !!
          snomlt = snomlt + snomlt_ros

          !! adjust for areal extent of snow cover
          if (hru(j)%sno_mm < hru(j)%sno%covmx) then
            rto_sno = hru(j)%sno_mm / hru(j)%sno%covmx
            snocov = rto_sno / (rto_sno + Exp(hru(j)%snocov1 - hru(j)%snocov2 * rto_sno))
          else
            snocov = 1.
          endif
          snomlt = snomlt * snocov
          if (snomlt < 0.) snomlt = 0.
          if (snomlt > hru(j)%sno_mm) snomlt = hru(j)%sno_mm
        endif !! end of melt calculation
        
        !!snowpack ripening and liquid water routing
        rain_mm = precip_eff
        total_liquid_input = snomlt + rain_mm
        
        hru(j)%sno_liq = hru(j)%sno_liq + total_liquid_input
        hru(j)%sno_mm = hru(j)%sno_mm - snomlt
          
        liq_capacity = hru(j)%sno_mm * liq_max_ratio
        
        if (hru(j)%sno_liq > liq_capacity) then
            !!snowpack collapses, releasing excess water
            precip_eff = hru(j)%sno_liq - liq_capacity
            hru(j)%sno_liq = liq_capacity
        else
            !!all liquid water retained inside the snowpack sponge
            precip_eff = 0.
        endif     
          
        if (time%step > 1) then
        w%ts(:) = w%ts(:) + snomlt / time%step
        end if
        if (precip_eff < 0.) precip_eff = 0.
        
      end if !!end of sno_mm > 0 block
      
      if (hru(j)%sno_mm < 1.e-6) then
          hru(j)%sno_tmp = 0. ! when there is no snowpack, then reset the temp of snowpack
          hru(j)%sno_liq = 0. !! Safeguard to clear liquid reservoir
      endif    
      return
      end subroutine sq_snom