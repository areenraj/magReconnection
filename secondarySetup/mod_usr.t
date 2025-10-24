module mod_usr
    use mod_mhd
    implicit none

    double precision :: q_e, unit_currentdensity
    double precision :: unit_electricfield
    double precision :: rhoConst, tempConst, pConst

    double precision :: parb, eta0, reta, heta, eta0_pre, reta_pre, heta_pre, tar, tacc

contains 

    subroutine usr_init()

        call set_coordinate_system("Cartesian_2.5D")

        unit_length        = 1.d9 ! cm
        unit_temperature   = 1.d6 ! K
        unit_numberdensity = 1.d9 ! cm^-3

        call usr_params_read(par_files)

        usr_init_one_grid       => initonegrid_usr
        usr_set_B0              => specialset_B0
        usr_set_J0              => specialset_J0
        usr_special_resistivity => special_eta
        usr_var_for_errest      => p_for_errest

        call mhd_activate()

        ! unit of current density
        unit_currentdensity=unit_magneticfield/unit_length/4.d0/dpi

        ! unit of electric field
        unit_electricfield=unit_magneticfield*unit_length/(unit_time*const_c)

        ! unit of charge
        q_e=unit_currentdensity/unit_numberdensity/unit_velocity

        if(mype==0) print*,'unit of charge',q_e

        ! dimensionless charge of electron
        q_e=1.60217653d-19/q_e

        if(mype==0) print*,'dimensionless e',q_e

        if (mype==0) then
        print *, 'rho', unit_density
        print *, 'B', unit_magneticfield
        print *, 'p', unit_pressure
        print *, 'v', unit_velocity
        print *, 't', unit_time
        print *, 'E', unit_electricfield
        print *, 'J', unit_currentdensity
        print *, 'eta', unit_electricfield/unit_currentdensity
        endif

    end subroutine usr_init

    subroutine usr_params_read(files)

        character(len=*), intent(in) :: files(:)
        integer                      :: n

        namelist /usr_list/ parb, &
            eta0_pre, reta_pre, heta_pre, tar, &
            eta0, reta, heta, tacc

    !    write(*,*) "begin usr_params_read",mype
        do n = 1, size(files)
        open(unitpar, file=trim(files(n)), status="old")
        read(unitpar, usr_list, end=111)
    111    close(unitpar)
    end do
    end subroutine usr_params_read

    subroutine initonegrid_usr(ixI^L,ixO^L,w,x)

        integer, intent(in) :: ixI^L, ixO^L
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: w(ixI^S,1:nw)

        rhoConst = 4.d5
        tempConst = 2.d6/unit_temperature
        pConst = tempConst * rhoConst

        w(ixO^S,rho_)= rhoConst
        w(ixO^S,mom(:)) = zero
        w(ixO^S,p_)= pConst
        w(ixO^S,mag(:))= zero

        call mhd_to_conserved(ixI^L,ixO^L,w,x)
    
    end subroutine initonegrid_usr

    subroutine p_for_errest(ixI^L,ixO^L,iflag,w,x,var)
        integer, intent(in)           :: ixI^L,ixO^L,iflag
        double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
        double precision, intent(out) :: var(ixI^S)

        call mhd_get_pthermal(w,x,ixI^L,ixO^L,var)
        
    end subroutine p_for_errest  

    subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
    ! Here add a time-independent background magnetic field
        integer, intent(in)           :: ixI^L,ixO^L
        double precision, intent(in)  :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: wB0(ixI^S,1:ndir)

        parb=1.d1

        wB0(ixO^S,1)=zero
        wB0(ixO^S,2)=-Busr*dtanh(parb*x(ixO^S,1))
        wB0(ixO^S,3)=Busr/dcosh(parb*x(ixO^S,1))

    end subroutine specialset_B0  

    subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
    ! Here add a time-independent background current density 
        integer, intent(in)           :: ixI^L,ixO^L
        double precision, intent(in)  :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)

        parb=1.d1

        wJ0(ixO^S,1)=zero
        wJ0(ixO^S,2)=parb*Busr*dtanh(parb*x(ixO^S,1))/dcosh(parb*x(ixO^S,1))
        wJ0(ixO^S,3)=-parb*Busr/dcosh(parb*x(ixO^S,1))**2

    end subroutine specialset_J0

    subroutine special_eta(w,ixI^L,ixO^L,idirmin,x,current,eta)
        
        integer, intent(in) :: ixI^L, ixO^L, idirmin
        double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)
        double precision :: current(ixI^S,7-2*ndir:3), eta(ixI^S)
        double precision :: rad(ixI^S)

        parb=1.d1
        eta0_pre=5.d-3
        reta_pre=0.048d0
        heta_pre=1.d0
        eta0=5.d-3
        reta=0.048d0
        heta=1.d0
        tar=1.67d0
        tacc=1.78d0

        if (global_time<tar) then
        rad(ixO^S)=dsqrt(x(ixO^S,1)**2+(x(ixO^S,2)-heta_pre)**2)
        where (rad(ixO^S) .lt. reta_pre)
            eta(ixO^S)=eta0_pre*(2.d0*(rad(ixO^S)/reta_pre)**3-3.d0*(rad(ixO^S)/reta_pre)**2+1.d0)
        elsewhere
            eta(ixO^S)=zero
        endwhere
        elseif (global_time<tacc) then
        rad(ixO^S)=dsqrt(x(ixO^S,1)**2+(x(ixO^S,2)-heta)**2)
        where (rad(ixO^S) .lt. reta)
            eta(ixO^S)=eta0*(2.d0*(rad(ixO^S)/reta)**3-3.d0*(rad(ixO^S)/reta)**2+1.d0)
        elsewhere
            eta(ixO^S)=zero
        endwhere
        else
            eta(ixO^S)=eta0
        end if 

    end subroutine special_eta   

end module mod_usr 