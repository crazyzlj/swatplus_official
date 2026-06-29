      subroutine stmp_solt
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine estimates daily average temperature at the bottom
!!    of each soil layer     

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    albday      |none          |albedo of ground for day
!!    tmp_an(:)   |deg C         |average annual air temperature
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Log, Max, Min

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
      use basin_module
      use climate_module
      use septic_data_module
      use hydrograph_module, only : ob,sp_ob,sp_ob1
      use hru_module, only : hru, iseptic, ihru, i_sep, iwgen, albday, isep 
      use soil_module
      use time_module
      use organic_mineral_mass_module
      use sd_channel_module
      
      implicit none

      integer :: j = 0           !none          |HRU number
      integer :: k = 0           !none          |counter
      integer :: iob = 0         !none          |object number
      integer :: iru = 0         !none          |routing unit object number
      integer :: icha = 0        !none          |sequential channel index
      real :: f = 0.             !none          |variable to hold intermediate calculation result
      real :: dp = 0.            !mm            |maximum damping depth
      real :: ww = 0.            !none          |variable to hold intermediate calculation
      real :: b = 0.             !none          |variable to hold intermediate calculation
      real :: wc = 0.            !none          |scaling factor for soil water impact on daily damping depth
      real :: dd = 0.            !mm            |damping depth for day
      real :: xx = 0.            !none          |variable to hold intermediate calculation
      real :: st0 = 0.           !MJ/m^2        |radiation hitting soil surface on day
      real :: tlag = 0.          !none          |lag coefficient for soil temperature
      real :: df = 0.            !none          |depth factor
      real :: zd = 0.            !none          |ratio of depth at center of layer to damping depth 
      real :: bcv = 0.           !none          |lagging factor for cover
      real :: tbare = 0.         !deg C         |temperature of bare soil surface
      real :: tcov = 0.          !deg C         |temperature of soil surface corrected for cover
      real :: tmp_srf = 0.       !deg C         |temperature of soil surface
      real :: snowpack_swe = 0.  !mm H2O        |solid snow plus retained liquid water
      real :: cover = 0.         !kg/ha         |soil cover
      real :: t_thaw = 0.        !deg C         |temperature at which frozen-state target is zero
      real :: t_froz = 0.        !deg C         |temperature at which frozen-state target is one
      real :: t_lyr2 = 0.        !deg C         |soil temperature used for frozen-state update
      real :: frz_target = 0.    !none          |instantaneous frozen-state target
      real :: frz_alpha = 0.

      j = ihru

      !tlag = 0.8
      tlag = bsn_prm%tlag

!! calculate damping depth

      !! calculate maximum damping depth
      !! SWAT manual equation 2.3.6
      f = 0.
      dp = 0.
      f = soil(j)%avbd / (soil(j)%avbd + 686. * Exp(-5.63 *       &       
              soil(j)%avbd))
      dp = 1000. + 2500. * f

      !! calculate scaling factor for soil water
      !! SWAT manual equation 2.3.7
      ww = 0.
      wc = 0.
      ww = .356 - .144 * soil(j)%avbd
      wc = soil(j)%sw / (ww * soil(j)%phys(soil(j)%nly)%d)

      !! calculate daily value for damping depth
      !! SWAT manual equation 2.3.8
      b = 0.
      f = 0.
      dd = 0.
      b = Log(500. / dp)
      f = Exp(b * ((1. - wc) / (1. + wc))**2)
      dd = f * dp

!! calculate lagging factor for soil cover impact on soil surface temp
!! SWAT manual equation 2.3.11
      cover = pl_mass(j)%ab_gr_com%m + pl_mass(j)%rsd_tot%m
      bcv = cover / (cover + Exp(7.563 - 1.297e-4 * cover))
      snowpack_swe = hru(j)%sno_mm + hru(j)%sno_liq
      if (snowpack_swe /= 0.) then
        if (snowpack_swe <= 120.) then
          xx = 0.
          xx = snowpack_swe / (snowpack_swe + Exp(6.055 - .3002 * snowpack_swe))
        else
          xx = 1.
        end if
        bcv = Max(xx,bcv)
      end if

!! calculate temperature at soil surface
      st0 = 0.
      tbare = 0.
      tcov = 0.
      tmp_srf = 0.
      !! SWAT manual equation 2.3.10
      st0 = (w%solrad * (1. - albday) - 14.) / 20.
      !! SWAT manual equation 2.3.9
      tbare = w%tave + 0.5 * (w%tmax - w%tmin) * st0
      !! SWAT manual equation 2.3.12
      tcov = bcv * soil(j)%phys(2)%tmp + (1. - bcv) * tbare

!!    taking average of bare soil and covered soil as in APEX
!!    previously using minimum causing soil temp to decrease
!!    in summer due to high biomass

      tmp_srf = 0.5 * (tbare + tcov)  ! following Jimmy"s code

!! calculate temperature for each layer on current day
      xx = 0.
      do k = 1, soil(j)%nly
        zd = 0.
        df = 0.
        zd = (xx + soil(j)%phys(k)%d) / 2.  ! calculate depth at center of layer
        zd = zd / dd                 ! SWAT manual equation 2.3.5
        !! SWAT manual equation 2.3.4
        df = zd / (zd + Exp(-.8669 - 2.0775 * zd))
        !! SWAT manual equation 2.3.3
        soil(j)%phys(k)%tmp = tlag * soil(j)%phys(k)%tmp + (1. - tlag) *       &
                      (df * (wgn_pms(iwgen)%tmp_an - tmp_srf) + tmp_srf)
        xx = soil(j)%phys(k)%d

        ! Temperature correction for Onsite Septic systems
        isep = iseptic(j)
        if (sep(isep)%opt /= 0 .and. time%yrc >= sep(isep)%yr .and. k >=       &
                                                          i_sep(j)) then
       if (soil(j)%phys(k)%tmp < 10.) then
           soil(j)%phys(k)%tmp = 10. - (10. - soil(j)%phys(k)%tmp) * 0.1
       end if     
      endif

      end do

      if (soil(j)%nly >= 2) then
        t_lyr2 = soil(j)%phys(2)%tmp
      else
        t_lyr2 = soil(j)%phys(1)%tmp
      endif

      if (bsn_cc%froz_soil == 0) then
          !! Original SWAT+ style: binary frozen condition using 0 deg C
          if (t_lyr2 <= 0.0) then
              soil(j)%frz_state = 1.0
          else
              soil(j)%frz_state = 0.0
          end if
      else
          !! Enhanced method: continuous frozen-soil state
          t_thaw = bsn_prm%frz_t_thaw
          t_froz = bsn_prm%frz_t_froz
          if (t_lyr2 >= t_thaw) then
            frz_target = 0.0
          else if (t_lyr2 <= t_froz) then
            frz_target = 1.0
          else
            frz_target = (t_thaw - t_lyr2) / (t_thaw - t_froz)
          end if
          frz_target = Max(0.0, Min(1.0, frz_target))

          if (frz_target > soil(j)%frz_state) then
              ! freezing or refreezing
              if (t_lyr2 <= -1.0) then
                  frz_alpha = bsn_prm%frz_alpha_fr_cold
              else
                  frz_alpha = bsn_prm%frz_alpha_fr_warm
              end if
          else
              ! thawing
              if (t_lyr2 >= 0.0) then
                  frz_alpha = bsn_prm%frz_alpha_th_warm
              else if (t_lyr2 >= -0.5) then
                  frz_alpha = bsn_prm%frz_alpha_th_cool
              else
                  frz_alpha = bsn_prm%frz_alpha_th_cold
              end if
          endif
          soil(j)%frz_state = soil(j)%frz_state + frz_alpha * &
                  (frz_target - soil(j)%frz_state)
          soil(j)%frz_state = Max(0.0, Min(1.0, soil(j)%frz_state))
      end if
      if (bsn_cc%icejam == 1) then
          !! Update frozen soil diagnostics for the downstream channel.
          iob = j + sp_ob1%hru - 1
          if (ob(iob)%ru_tot > 0) then
              iru = ob(iob)%ru(1)                     ! lsu number
              iru = iru + sp_ob1%ru - 1
              if (ob(iru)%src_tot > 0) then
                  if (ob(iru)%obtyp_out(1) == 'sdc') then
                      icha = ob(iru)%obj_out(1)          ! channel object index
                      icha = icha - sp_ob1%chandeg + 1  ! sequential channel index
                      if (icha > 0 .and. icha <= size(sd_ch)) then
                          sd_ch(icha)%frz_surf_avg = sd_ch(icha)%frz_surf_avg + &
                                  soil(j)%frz_state ** bsn_prm%frz_surf_exp * hru(j)%area_ha
                          if (soil(j)%frz_state >= 0.5) then
                              sd_ch(icha)%frz_area_frac = sd_ch(icha)%frz_area_frac + hru(j)%area_ha
                          endif
                      endif
                  endif
              endif
          endif
      end if
      
      return
      end subroutine stmp_solt