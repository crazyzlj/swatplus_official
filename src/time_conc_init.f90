      subroutine time_conc_init 
    
      use ru_module
      use hru_module, only : brt, hru, ihru, t_ov, tconc
      use hydrograph_module, only : sp_ob, ru_def, ru_elem, sp_ob1, ob
      use topography_data_module
      use time_module
      use basin_module
      
      implicit none 
      
      !!compute time of concentration for routing units
      do iru = 1, sp_ob%ru
          call ru_tc_upd(iru)
      enddo
      
      !!compute time of concentration (sum of overland and channel times) for hru
      do ihru = 1, sp_ob%hru
          call compute_hru_routing(ihru, hru(ihru)%luse%ovn, bsn_prm%surlag)
      enddo    
      return
      end subroutine time_conc_init
      
      subroutine ru_tc_upd(j) 
      use ru_module
      use hru_module, only : hru, ihru
      use hydrograph_module, only : sp_ob, ru_def, ru_elem, sp_ob1, ob
      use soil_module, only : soil
    
      implicit none 
      
      integer, intent(in) :: j
      
      integer :: ii = 0            !none         |counter
      integer :: ielem = 0         !none         |counter
      integer :: iob = 0           !             | 
      real :: current_ovn_hru = 0. !none        |manning's roughness coefficient of hru
    
     ! compute weighted Mannings n for each subbasin
      !do iru = 1, sp_ob%ru
        iru = j
        ru_n(iru) = 0.
        do ii = 1, ru_def(iru)%num_tot
          ielem = ru_def(iru)%num(ii)
          if (ru_elem(ielem)%obtyp == "hru") then
            ihru = ru_elem(ielem)%obtypno 
            current_ovn_hru = hru(ihru)%luse%ovn * (1.0 - soil(ihru)%frz_state) + (0.01 * soil(ihru)%frz_state)
            ru_n(iru) = ru_n(iru) + current_ovn_hru * hru(ihru)%km
          else
            ru_n(iru) = 0.1
          end if
        end do
        iob = sp_ob1%ru + iru - 1
        ru(iru)%da_km2 = ob(iob)%area_ha / 100.
        ru_n(iru) = ru_n(iru) / ru(iru)%da_km2
        
        call compute_ru_routing(iru, ru_n(iru))
      !end do
      end subroutine ru_tc_upd
            
      subroutine compute_ru_routing(j, run)
        use ru_module
        use topography_data_module
        
        implicit none
        
        integer, intent(in) :: j
        real, intent(in) :: run
        
        integer :: ith = 0           !             |
        integer :: ifld = 0          !             |
        real :: tov = 0.             !             |
        real :: ch_slope = 0.        !             |
        real :: ch_n = 0.            !             |
        real :: ch_l = 0.            !             | 
        real :: t_ch = 0.            !hr           |time for flow entering the farthest upstream 
                                     !             |channel to reach the subbasin outlet
        ith = ru(j)%dbs%toposub_db
        !if (ith > 0 .and. ichan > 0) then                  
        ! compute tc for the subbasin
          tov = .0556 * (topo_db(ith)%slope_len * run) ** .6 /     &
                                              (topo_db(ith)%slope + .001) ** .3
          ch_slope = .5 * (topo_db(ith)%slope + .001)
          ch_n = run
          ch_l = ru(j)%field%length / 1000.
          t_ch = .62 * ch_l * ch_n**.75 / (ru(j)%da_km2**.125 * ch_slope**.375)
          ru_tc(j) = tov + t_ch
        !end if                                             
      end subroutine compute_ru_routing
      
    
      subroutine time_conc_upd(j) 
    
      use hru_module, only : hru
      use basin_module, only: bsn_prm
      use soil_module, only : soil
      
      implicit none 
      
      integer, intent(in) :: j    !none   |HRU number
      real :: current_ovn         !none   |manning's roughness coefficient   
      real :: current_surlag      !days   |lag
      
      current_ovn = hru(j)%luse%ovn * (1.0 - soil(j)%frz_state) + (0.01 * soil(j)%frz_state)
      current_surlag = bsn_prm%surlag * (1.0 - soil(j)%frz_state) + (24.0 * soil(j)%frz_state)

      call compute_hru_routing(j, current_ovn, current_surlag)
      
      return
      end subroutine time_conc_upd
    
      subroutine compute_hru_routing(ihru, input_ovn, input_surlag)
      
      use hru_module, only : hru, t_ov, tconc, brt
      use topography_data_module, only : topo_db
      use time_module, only : time
      
      implicit none
      
      integer, intent(in) :: ihru         !none   |HRU number
      real, intent(in) :: input_ovn       !none   |manning's roughness coefficient
      real, intent(in) :: input_surlag    !days   |lag
      
      integer :: ith = 0
      real :: ch_slope = 0.
      real :: ch_l = 0.
      real :: t_ch = 0.
      
        ith = hru(ihru)%dbs%topo
        t_ov(ihru) = .0556 * (hru(ihru)%topo%slope_len *                    &
           input_ovn) ** .6 / (hru(ihru)%topo%slope + .0001) ** .3
        ch_slope = .5 * topo_db(ith)%slope
        !ch_n = hru(ihru)%luse%ovn
        !! assume length to width (l/w) ratio of 2 --> A=l*w - A=l*l/2 - l=sqrt(A/2)
        !! assume channel begins at 1/2 of distance
        ch_l = 0.5 * sqrt(hru(ihru)%area_ha / 2.)
        !ch_l = hru(ihru)%field%length / 1000.
        t_ch = .31 * ch_l * input_ovn**.75 / (hru(ihru)%km**.125 * (ch_slope + .001)**.375)
        tconc(ihru) = t_ov(ihru) + t_ch
        !! compute fraction of surface runoff that is reaching the main channel
        if (time%step > 1) then
          brt(ihru) = 1.-Exp(-input_surlag / (tconc(ihru) /               &
              (time%dtm / 60.)))    !! urban modeling by J.Jeong
        else
          brt(ihru) = 1. - Exp(-input_surlag / tconc(ihru))
        endif
      return
      end subroutine compute_hru_routing
      
     