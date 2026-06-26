      subroutine dr_ru
    
      use hydrograph_module
      use hru_lte_module
      use ru_module
      use hru_module, only : ihru, tconc
      
      implicit none

      integer :: ii = 0             !none          |counter
      
     ! compute delivery ratio for each hru in the sub
      do ii = 1, sp_ob%ru
          call dr_ru_upd(ii)
      end do

      return
      end subroutine dr_ru
          
      subroutine dr_ru_upd(jru)
      
      use hydrograph_module
      use ru_module
      
      implicit none
      
      integer, intent(in) :: jru
      integer :: ii = 0
      integer :: ielem = 0
      
      do ii = 1, ru_def(jru)%num_tot
        ielem = ru_def(jru)%num(ii)
        call dr_elem_upd(jru, ielem)
      end do  
      return
      end subroutine dr_ru_upd    
      
      subroutine dr_elem_upd(jru, ielem)

      use hydrograph_module
      use hru_lte_module
      use ru_module
      use hru_module, only : ihru, tconc

      implicit none

      integer, intent(in) :: jru
      integer, intent(in) :: ielem
      real :: rto = 0.

      if (ru_elem(ielem)%dr_name == "calc" .or. ru_elem(ielem)%dr_name == "0") then
        select case (ru_elem(ielem)%obtyp)
        case ("hru")
            ihru = ru_elem(ielem)%obtypno
            if (ru_tc(jru) > 1.e-6) then
                rto = tconc(ihru) / ru_tc(jru)
            else
                rto = 1.
            endif
        case ("hlt")
            ihru = ru_elem(ielem)%obtypno
            if (ru_tc(jru) > 1.e-6) then
                rto = (hlt_db(ihru)%tc / 3600.) / ru_tc(jru)
            else
                rto = 1.
            endif
        case ("sdc")
            rto = 1.
        case ("ru")
            rto = 1.
        case default
            rto = 1.
        end select
            
        rto = amin1(1.0, rto ** .5)
        ru_elem(ielem)%dr = rto .add. hz
        ru_elem(ielem)%dr%flo = 1.
      end if
        
      if (ru_elem(ielem)%dr_name == "full" .or. ru_elem(ielem)%dr_name == "1") then
        ru_elem(ielem)%dr = 1. .add. hz
      end if
      
      return
      end subroutine dr_elem_upd   