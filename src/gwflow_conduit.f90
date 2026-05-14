      subroutine gwflow_conduit(chan_id) !ljzhu 05/11/2026

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates the water discharge volume from groundwater grid cells to the connected channel
!!    (exchange volumes are used in gwflow_simulate, in groundwater balance equations)

      use gwflow_module
      use hydrograph_module, only : ch_stor,ch_out_d, hz
      use constituent_mass_module

      implicit none

      integer, intent (in) :: chan_id    !       |channel number
      integer :: k = 0                   !       |counter for cells connected to the channel
      integer :: cell_id = 0             !       |id of cell connected to the channel
      real :: stor_volume = 0.           !m3     |
      real :: cdut_elev = 0.             !m      |elevation of conduit
      real :: head_diff = 0.             !m      |head difference between groundwater head and conduit
      real :: Q_depth = 0.               !mm/day |
      real :: Q = 0.                     !m3/day |  
      real :: Qup = 0.                   !m3/day |
      real :: Qchan = 0.                 !m3/day |     
      real :: Qleak = 0.                 !m3/day |      
      real :: excess = 0.


      !only proceed if conduit is active
      if (gw_conduit_flag == 1) then
        gw_conduit_info(chan_id)%output = hz ! reset 

        !loop through the cells connected to the channel
        do k=1,gw_conduit_info(chan_id)%ncon

          !cell in connection with the channel
          cell_id = gw_conduit_info(chan_id)%cells(k)

          !only proceed if the cell is active
          if(gw_state(cell_id)%stat == 1) then
            !excess = max(0., gw_cdut_stor(cell_id) - gw_cdut_smin(cell_id))
            !Qup = gw_cdut_k(cell_id) * excess ** gw_cdut_exp(cdut_id)
            !Qup = min(Qup, gw_cdut_qmax(cell_id))
            !Qup = min(Qup, gw_cdut_stor(cell_id))
            !Qchan = (1. - gw_cdut_leak(cell_id)) * Qup
            excess = max(0., gw_cdut_stor(cell_id) - 10.)
            Qup = 0.8 * excess
            !Qup = min(Qup, gw_cdut_qmax(cell_id))
            Qup = min(Qup, gw_cdut_stor(cell_id))
            Qchan = (1. - 0.1) * Qup  
            Qleak = Qup - Qchan
            
            gw_cdut_stor(cell_id) = gw_cdut_stor(cell_id) - Qup
  
            gw_conduit_info(chan_id)%output%flo = gw_conduit_info(chan_id)%output%flo + Qchan
            
            gw_hyd_ss(cell_id)%cdut = gw_hyd_ss(cell_id)%cdut + Qleak !entering aquifer
            gw_hyd_ss_yr(cell_id)%cdut = gw_hyd_ss_yr(cell_id)%cdut + Qleak !entering aquifer - store for annual water
            gw_hyd_ss_mo(cell_id)%cdut = gw_hyd_ss_mo(cell_id)%cdut + Qleak !entering aquifer - store for monthly water
              
            !!get head difference between groundwater head and conduit elevation (=channel bed elevation + gw_cdut_depth, calculated in gwflow_read)
            !head_diff = gw_state(cell_id)%head - (gw_state(cell_id)%elev - gw_cdut_depth(cell_id)) !m
            ! 
            !!only perform calculation if water table is above the conduit elevation
            !if(head_diff > 0.) then
            !  Q_depth = gw_cdut_conddepth(cell_id) * head_diff ** gw_cdut_exp(cell_id) !mm/day
            !  Q_depth = min(Q_depth, gw_cdut_qmax(cell_id))
            !  Q = Q_depth * gw_state(cell_id)%area / 1000. !m3/day
            !  Q = min(Q, gw_state(cell_id)%stor) !check for available groundwater in the cell - can only remove what is there
            ! 
            !  gw_state(cell_id)%stor = gw_state(cell_id)%stor - Q !update available groundwater in the cell
            !  gw_hyd_ss(cell_id)%cdut = gw_hyd_ss(cell_id)%cdut + Q * (-1) !leaving aquifer
            !  gw_hyd_ss_yr(cell_id)%cdut = gw_hyd_ss_yr(cell_id)%cdut + (Q*(-1)) !leaving aquifer - store for annual water
            !  gw_hyd_ss_mo(cell_id)%cdut = gw_hyd_ss_mo(cell_id)%cdut + (Q*(-1)) !leaving aquifer - store for monthly water
            !
            !  !add water to channel
            !  !ch_stor(chan_id)%flo = ch_stor(chan_id)%flo + Q !do not add to ch_stor, since it has little impact on the flow out rate.
            !  gw_conduit_info(chan_id)%output%flo = gw_conduit_info(chan_id)%output%flo + Q
            !
            !endif !check if groundwater head above conduit
          endif !check if cell is active
        enddo !go to next channel
        
        !if (cell_id == 3606) then
        !    write (9003,*) "after conduit, Q:", Q
        !endif  
        if (chan_id == 129) then
          stor_volume = gw_cdut_stor(2824) + gw_cdut_stor(3489) + gw_cdut_stor(3606) + gw_cdut_stor(3730) + gw_cdut_stor(3964)
          write (9003,*) "after conduit, cdut_stor:", stor_volume, "to channel:", gw_conduit_info(chan_id)%output%flo
        endif
      endif !check if conduit is active

      return
    end subroutine gwflow_conduit
