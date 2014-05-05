program main
    use omp_lib
    implicit none

    real, parameter :: packingFraction=0.2 ! freezing for HS at 0.494
    integer, parameter :: nst=10000 ! nst is the number of Monte Carlo steps
    integer, parameter :: np=10000 ! np == number of particles
    real, parameter :: dia=1. ! particles have diameter = 1.
    logical, parameter :: log_printpositionsXYZ = .false. ! print every positions at every MC step
    integer, parameter :: dimensionality=3
    real, allocatable :: x(:,:), xold(:,:), distSQarray(:), histo(:), randomRealArray(:,:)
    ! x(i,j) is the position of particle j along dir i, in lab frame
    real, parameter :: pi=acos(-1.)
    integer :: i, chosenParticle, nref,j
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
        stop "dimensionaly must be 3"
    end select

    nbinhistogram=int(lx/2.)*100 ! histograms will have a resolution of 1/100th of HS dia

    allocate( x              (dimensionality,np) ,source=0.)
    allocate( xold           (dimensionality,np) ,source=0.)
    allocate( randomRealArray(dimensionality,np) ,source=0.)

    call generateStartingConfigurationWithLowEnergy
    en0 = sysEn(x) ! energy of the initial configuration

    nref= 0
    randomMoveScaler = dia/10000. ! magic number
    do i =1,nst
        if (log_printpositionsXYZ) call printpositionsXYZ
        if (modulo(i,nst/10)==0) print*,"step ",i,"/",nst
!~         if (modulo(i,freqadjustRandomMoveScaler)==0) call adjustRandomMoveScaler !DOES NOT WORK
        call random_number(randomRealArray)
        randomRealArray = (randomRealArray-0.5)*randomMoveScaler ! in [0.,1.[ => ()-0.5 in [-0.5,0.5[ => ()*lx/100. in [-0.05,0.05[   OK
        xold = x ! we backup x of the previous case in case of step rejection
        x = x + randomRealArray ! update positions
        where (x>=lx) x=x-lx ! for cubic supercells only
        where (x<0.) x=x+lx
        en1 = sysEn(x) ! new energy
        dE= en1-en0 ! difference between new and old energies
        if ( dE>0.) then
            call random_number(dice)
            if ( dice >= exp(-dE) ) then
                nref = nref+1
                x = xold
                en1 = en0
            end if
        end if
        distSQarray = sum((x-lx/2.)**2,1) ! table of distances to supercell center
        call psl_histogram(rmin=0.,rmax=lx/2.,nbin=nbinhistogram,array1d=distSQarray,histo=histo)
    end do
print*,"acceptance ratio =",nst-nref,"/",nst,"=",real(nst-nref)/real(nst)
!~     call psl_normhistogram(histo)
    call psl_printhistogram(rmin=0.,rmax=lx/2.,nbin=nbinhistogram,histo=histo)

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
        real :: distsq
        integer :: i, j
        sysEn = 0.
        
        !solvent-solvent
        do i=1,np-1
            do j=i+1,np
                distsq =  sum( (x(:,i)-x(:,j))**2 )
                if ( distsq < dia**2 ) then
                    sysEn = huge(1.)
                    return
                end if
            end do
        end do

        !solvent-solute
        do i=1,np
            distsq = sum( (x(:,i)-lx*[0.5,0.5,0.5])**2 ) ! distance between solute at supercell center and solvent
            if ( distsq < dia**2 ) then
                sysEn = huge(1.)
                return
            end if
        end do

    end function sysEn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine generateStartingConfigurationWithLowEnergy
        implicit none
        integer :: i,n,j,k
        integer :: maximumNumberOfTries=1000
        logical :: doitagain
        real :: t0,t1,r(3)
        call cpu_time(t0)
        print*,"=> Random supercell initialization ..."
        i=0
        do while (sysEn(x)>epsilon(1.) .and. i<=maximumNumberOfTries)
            i=i+1
            call random_number(x)
            x=x*lx
        end do
        print*,"=> Random supercell initialization by method 1 failed."

        x=0.
        n=1
        do while (n<np+1)
            doitagain=.false.
            call random_number(x(:,n))
            x(:,n)=x(:,n)*lx
            if (sum( (x(:,n)-lx/2.)**2 )<dia**2) cycle
            do concurrent(i=1:np, i/=n)
                if (sum((x(:,n)-x(:,i))**2)<dia**2) doitagain=.true.
            end do
            if (.not. doitagain) n=n+1
        end do

        print*,"=> Random supercell initialization by method 2 done."

        if (sysen(x)>epsilon(1.)) stop "System energy not 0. Abord"
        if (any(x>=lx)) stop "Some position is outside [0,lx["

        call cpu_time(t1)
        print*,"=> Random supercell initialization took ",t1-t0,"sec"
    end subroutine generateStartingConfigurationWithLowEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine adjustRandomMoveScaler
        implicit none
        real, parameter :: targetAcceptanceRatio=0.3 ! papier K. Binder
        real :: acceptanceRatio
        STOP "DOES NOT WORK"
        acceptanceRatio = 1.-real(nref)/real(freqadjustRandomMoveScaler) ! in [0,1]
print*,"acceptanceRatio before=",acceptanceRatio
        if (acceptanceRatio<=epsilon(1.)) then
            randomMoveScaler = dia/10.
        else
            randomMoveScaler = randomMoveScaler * acceptanceRatio/targetAcceptanceRatio/100.
!~             randomMoveScaler = randomMoveScaler * min(1.1,max(0.9 , acceptanceRatio/targetAcceptanceRatio))
        end if
print*,"acceptanceRatio after=",acceptanceRatio
        nref = 0
    end subroutine adjustRandomMoveScaler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine printpositionsXYZ
        implicit none
        write(88,*)np
        write(88,*)
        do j=1,np
            write(88,*)"H",x(:,j)
        end do
    end subroutine printpositionsXYZ

end program main
