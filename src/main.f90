program main
    use omp_lib
    implicit none

    real, parameter :: packingFraction=0.25 ! freezing for HS at 0.494
    real, parameter :: dia=1. ! particles have diameter = 1.
    integer, parameter :: nst=10000 ! nst is the number of Monte Carlo steps
    integer, parameter :: np=5 ! np == number of particles
    real, allocatable :: x(:,:), histo(:), randomRealArray(:), randomMoveArray(:)
    ! x(i,j) is the position of particle j along dir i, in lab frame
    real, parameter :: pi=acos(-1.)
    integer, parameter :: dimensionality=3
    integer :: i, chosenParticle, nref
    real :: en0, en1 ! total energy of the system at step i and i+1. The acceptance depends upon en1-en0
    real :: randomReal, randomMove, dice, dE, randomMoveScaler
    real :: lx ! lx is the length of the well in units of hard sphere diameters
    integer :: nbinhistogram
    integer, parameter :: freqadjustRandomMoveScaler=10000

    select case (dimensionality)
    case (3)
        block
            real :: vsphere
            vsphere = pi/6.*dia**3
            lx=(real(np)*vsphere/packingFraction)**(1./3.)
        end block
        print*,"************************************"
        print*,"packing fraction = ",packingFraction
        print*,"density (N/V) = ",real(np)/lx**3
        print*,"box length = ",lx
        print*,"number of particles = ",np
        print*,"************************************"
    case default
        stop "dimensionaly not 3"
    end select

    nbinhistogram=int(lx)*100

    allocate( x(dimensionality,np) ,source=0.)
    allocate( randomMoveArray(np) ,source=0.)
    allocate( randomRealArray(np) ,source=0.)

    call generateStartingConfigurationWithLowEnergy

    en0 = sysEn(x)
    nref= 0
    randomMoveScaler = dia/10. ! magic number
    do i =1,nst
        if (modulo(i,nst/10)==0) print*,"step ",i,"/",nst
        if (modulo(i,freqadjustRandomMoveScaler)==0) call adjustRandomMoveScaler
        call random_number(randomRealArray)
        randomMoveArray = (randomRealArray-0.5)*randomMoveScaler ! in [0.,1.[ => ()-0.5 in [-0.5,0.5[ => ()*lx/100. in [-0.05,0.05[   OK
        x(1,:) = x(1,:) + randomMoveArray ! update positions
        where (x(1,:)>=lx) x(1,:)=x(1,:)-lx
        where (x(1,:)<0) x(1,:)=x(1,:)+lx
        en1 = sysEn(x)
        dE= en1-en0
        if ( dE>0.) then
            call random_number(dice)
            if ( dice >= exp(-dE) ) then
                x(1,:) = x(1,:) - randomMoveArray
                nref = nref+1
                en1 = en0
            end if
        end if
!~         call psl_histogram(rmin=0.,rmax=lx,nbin=nbinhistogram,array1d=x(1,:),histo=histo)
    end do

!~     call psl_normhistogram(histo)
!~     call psl_printhistogram(rmin=0.,rmax=lx,nbin=nbinhistogram,histo=histo)

    contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine psl_printhistogram(rmin,rmax,nbin,histo)
        implicit none
        real, intent(in) :: rmin, rmax, histo(:)
        integer, intent(in) :: nbin
        integer :: i
        real :: bornmin, bornmax, dr
        if (ubound(histo,1)/=nbin) stop "upper bound of histo is not nbin"
        if (lbound(histo,1)/=1) stop "lower bound of histo should be 1"
        dr = (rmax-rmin)/real(nbin)
        do i=1,nbin
            bornmin =real(i-1)*dr
            bornmax =real(i)*dr
            write(77,*) bornmin, histo(i)
            write(77,*) bornmax, histo(i)
        end do
    end subroutine psl_printhistogram

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine psl_histogram (rmin,rmax,nbin,array1d,histo)
        implicit none
        real, intent(in) :: rmin, rmax, array1d(:)
        real, allocatable, intent(inout) :: histo(:)
        integer, intent(in) :: nbin
        integer :: i
        real :: bornmin, bornmax, dr, toadd
        if (.not. allocated(histo)) allocate(histo(nbin) ,source=0.)
        dr = (rmax-rmin)/real(nbin)
        do i=1,nbin
            bornmin =real(i-1)*dr
            bornmax =real(i)*dr
            toadd = real(count( array1d>=bornmin .and. array1d<bornmax ))  ! in [bornmin,bornmax[
            histo(i) = histo(i) + toadd
        end do
    end subroutine psl_histogram

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine psl_normhistogram (histo)
        implicit none
        real, intent(inout) :: histo(:)
        histo = histo/maxval(histo)
    end subroutine psl_normhistogram

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! hard sphere repulsion between particles of radius dia
    pure real function sysEn(x)
        implicit none
        real, intent(in) :: x(:,:)
        real :: dist
        integer :: i, j
        sysEn = 0.
        do i=1,np-1
            do j=i+1,np
                dist = abs(x(1,i)-x(1,j))
                if ( dist < dia ) then
                    sysEn = huge(1.)
                    return
                end if
            end do
        end do
    end function sysEn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine generateStartingConfigurationWithLowEnergy
        implicit none
        integer :: i
        print*,"=> supercell initialization"
        i=0
        do while (sysEn(x)>epsilon(1.))
            i=i+1
            if( i==1000 ) stop "I tried 1000 times to generate initial config. Now stop"
            call random_number(x(1,:))
            x=x*lx
        end do
        print*,"=> supercell initialization OK"
    end subroutine generateStartingConfigurationWithLowEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine adjustRandomMoveScaler
        implicit none
        real, parameter :: targetAcceptanceRatio=0.3 ! papier K. Binder
        real :: acceptanceRatio
        acceptanceRatio = 1.-real(nref)/real(freqadjustRandomMoveScaler)
        randomMoveScaler = randomMoveScaler * min(1.1,max(0.9 , acceptanceRatio/targetAcceptanceRatio))
        nref = 0
    end subroutine adjustRandomMoveScaler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program main
