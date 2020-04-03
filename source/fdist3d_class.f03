! fdist3d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module fdist3d_class

      use perrors_class
      use parallel_pipe_class
      use spect3d_class
      use ufield3d_class
      use part3d_lib
      use input_class
         
      implicit none

      private

      public :: fdist3d, fdist3d_000, fdist3d_001, fdist3d_002, fdist3d_100
      public :: fdist3d_003


      type, abstract :: fdist3d

         private

         class(spect3d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
!
! ndprof = profile type
         integer :: npf, npmax
         logical :: evol = .true.
                         
         contains
         
         generic :: new => init_fdist3d         
         generic :: del => end_fdist3d
         generic :: dist => dist3d
         procedure(ab_init_fdist3d), deferred, private :: init_fdist3d !输入文件写参数
         procedure, private :: end_fdist3d
         procedure(ab_dist3d), deferred, private :: dist3d !粒子初始化
         procedure :: getnpf, getnpmax, getevol
                  
      end type 

      abstract interface
!
      subroutine ab_dist3d(this,part3d,npp,fd)
         import fdist3d
         import ufield3d
         implicit none
         class(fdist3d), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part3d
         integer, intent(inout) :: npp
         class(ufield3d), intent(in), pointer :: fd         
      end subroutine ab_dist3d
!
      subroutine ab_init_fdist3d(this,input,i)
         import fdist3d
         import input_json
         implicit none
         class(fdist3d), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i
      end subroutine ab_init_fdist3d
!
      end interface
!
      type, extends(fdist3d) :: fdist3d_000
! Tri Gaussian profile (the same particle charge)
         private

         integer :: npx, npy, npz !粒子数
         real :: qm, sigx, sigy, sigz !荷质比
         real :: bcx, bcy, bcz, sigvx, sigvy, sigvz !beam center,vz 能散
         real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np !重心 不稳定性 重心能量 
         logical :: quiet !减少数值噪声

         contains
         procedure, private :: init_fdist3d => init_fdist3d_000
         procedure, private :: dist3d => dist3d_000 !赋值

      end type fdist3d_000
!
      type, extends(fdist3d) :: fdist3d_001
! Bi Gaussian and piecewise in z (the same particle charge)
         private

         integer :: npx, npy, npz
         real :: qm, sigx, sigy
         real :: bcx, bcy, bcz, sigvx, sigvy, sigvz
         real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np
         real, dimension(:), allocatable :: fz, z !分段分布
         logical :: quiet

         contains
         procedure, private :: init_fdist3d => init_fdist3d_001
         procedure, private :: dist3d => dist3d_001

      end type fdist3d_001
!
      type, extends(fdist3d) :: fdist3d_002
! Twiss parameter and Gaussian in z (the same particle charge)
         private

         integer :: npx, npy, npz
         real :: qm, sigvz, sigz
         real :: bcx, bcy, bcz, alphax, alphay, betax, betay
         real :: emitx, emity
         real :: gamma,np
         real, dimension(:), allocatable :: fz, z
         logical :: quiet

         contains
         procedure, private :: init_fdist3d => init_fdist3d_002
         procedure, private :: dist3d => dist3d_002

      end type fdist3d_002
!
      type, extends(fdist3d) :: fdist3d_003
! Read external particle data files
         private

         integer :: npt
         real :: bcx, bcy, bcz, dx, dy, dz, cwp, qm
         character(len=:), allocatable :: file !文件名

         contains
         procedure, private :: init_fdist3d => init_fdist3d_003
         procedure, private :: dist3d => dist3d_003

      end type fdist3d_003
!
      type, extends(fdist3d) :: fdist3d_100
! Ring profile
         private

         integer :: xppc, yppc, zppc
         real :: qm, bcx, bcy, bcz
         real, dimension(3) :: lb, ub
         real :: r1, r2
         real :: gamma, np

         contains
         procedure, private :: init_fdist3d => init_fdist3d_100
         procedure, private :: dist3d => dist3d_100

      end type fdist3d_100

      character(len=10), save :: class = 'fdist3d:'
      character(len=128), save :: erstr
      
      contains
!
      function getnpf(this)

         implicit none

         class(fdist3d), intent(in) :: this
         integer :: getnpf
         
         getnpf = this%npf

      end function getnpf
!      
      function getnpmax(this)

         implicit none

         class(fdist3d), intent(in) :: this
         integer :: getnpmax
         
         getnpmax = this%npmax

      end function getnpmax
!      
      function getevol(this)

         implicit none

         class(fdist3d), intent(in) :: this
         logical :: getevol
         
         getevol = this%evol

      end function getevol
!      
      subroutine end_fdist3d(this)
          
         implicit none
         
         class(fdist3d), intent(inout) :: this
         character(len=18), save :: sname = 'end_fdist3d:'

         call this%err%werrfl2(class//sname//' started')
         call this%err%werrfl2(class//sname//' ended')
                  
      end subroutine end_fdist3d
!
      subroutine init_fdist3d_000(this,input,i) !000读参数
      
         implicit none
         
         class(fdist3d_000), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i        
! local data
         integer :: npf,npx,npy,npz,npmax
         real :: qm,sigx,sigy,sigz,bcx,bcy,bcz,sigvx,sigvy,sigvz
         real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np
         logical :: quiet, evol
         real :: min, max, cwp, n0
         real :: alx, aly, alz, dx, dy, dz
         integer :: indx, indy, indz         
         character(len=20) :: sn,s1
         character(len=18), save :: sname = 'init_fdist3d_000:'
         
         this%sp => input%sp
         this%err => input%err
         this%p => input%pp

         call this%err%werrfl2(class//sname//' started')

         write (sn,'(I3.3)') i
         s1 = 'beam('//trim(sn)//')' !！几个beam

         call input%get('simulation.n0',n0)
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get('simulation.indz',indz)

         cwp=5.32150254*1e9/sqrt(n0)
         call input%get('simulation.box.x(1)',min)
         call input%get('simulation.box.x(2)',max)
         call input%get(trim(s1)//'.center(1)',bcx)
         bcx = bcx - min
         alx = (max-min) 
         dx=alx/real(2**indx)
         call input%get('simulation.box.y(1)',min)
         call input%get('simulation.box.y(2)',max)
         call input%get(trim(s1)//'.center(2)',bcy)
         bcy = bcy -min
         aly = (max-min) 
         dy=aly/real(2**indy)
         call input%get('simulation.box.z(1)',min)
         call input%get('simulation.box.z(2)',max)
         call input%get(trim(s1)//'.center(3)',bcz)
         bcz = bcz -min
         alz = (max-min) 
         dz=alz/real(2**indz)

         call input%get(trim(s1)//'.profile',npf) !从输入文件读
         call input%get(trim(s1)//'.np(1)',npx)
         call input%get(trim(s1)//'.np(2)',npy)
         call input%get(trim(s1)//'.np(3)',npz)
         call input%get(trim(s1)//'.q',qm)
         call input%get(trim(s1)//'.sigma(1)',sigx)
         call input%get(trim(s1)//'.sigma(2)',sigy)
         call input%get(trim(s1)//'.sigma(3)',sigz)
         call input%get(trim(s1)//'.sigma_v(1)',sigvx)
         call input%get(trim(s1)//'.sigma_v(2)',sigvy)
         call input%get(trim(s1)//'.sigma_v(3)',sigvz)
         call input%get(trim(s1)//'.centroid_x(1)',cx1)
         call input%get(trim(s1)//'.centroid_x(2)',cx2)
         call input%get(trim(s1)//'.centroid_x(3)',cx3)
         call input%get(trim(s1)//'.centroid_y(1)',cy1)
         call input%get(trim(s1)//'.centroid_y(2)',cy2)
         call input%get(trim(s1)//'.centroid_y(3)',cy3)
         call input%get(trim(s1)//'.quiet_start',quiet)
         call input%get(trim(s1)//'.gamma',gamma)
         call input%get(trim(s1)//'.peak_density',np)
         call input%get(trim(s1)//'.npmax',npmax)
         call input%get(trim(s1)//'.evolution',evol)



         this%npf = npf !单位变换
         this%npx = npx
         this%npy = npy
         this%npz = npz
         this%npmax = npmax
         qm = qm/abs(qm)*abs(np)*(2*3.1415926535897932)**1.5*sigx*sigy*sigz
         qm = qm/dx/dy/dz
         qm = qm/npx
         qm = qm/npy
         qm = qm/npz
         this%qm = qm
         this%bcx = bcx/dx
         this%bcy = bcy/dy
         this%bcz = bcz/dz
         this%sigx = sigx/dx
         this%sigy = sigy/dy
         this%sigz = sigz/dz
         this%sigvx = sigvx
         this%sigvy = sigvy
         this%sigvz = sigvz
         this%cx1 = cx1*dz*dz/dx
         this%cx2 = cx2*dz/dx
         this%cx3 = cx3/dx
         this%cy1 = cy1*dz*dz/dy
         this%cy2 = cy2*dz/dy
         this%cy3 = cy3/dy
         this%gamma = gamma
         this%np = np
         this%quiet = quiet
         this%evol = evol

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist3d_000
!
      subroutine dist3d_000(this,part3d,npp,fd)
      
         implicit none
         
         class(fdist3d_000), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part3d
         integer, intent(inout) :: npp
         class(ufield3d), intent(in), pointer :: fd
! local data1

! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
         real, dimension(:,:), pointer :: pt => null()
         integer :: npx, npy, npz, nx, ny, nz, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real :: sigx, sigy, sigz, x0, y0, z0
         real, dimension(3) :: cx, cy
         real, dimension(4) :: edges
         integer, dimension(2) :: noff
         integer :: nps=1
         logical :: lquiet = .false.
         integer :: idimp, npmax, ierr = 0
         character(len=18), save :: sname = 'dist3d_000:'

         call this%err%werrfl2(class//sname//' started')
         
         npx = this%npx; npy = this%npy; npz = this%npz
         nx = fd%getnd1(); ny = fd%getnd2(); nz = fd%getnd3()
         ipbc = this%sp%getpsolver()
         pt => part3d
         vtx = this%sigvx; vty = this%sigvy; vtz = this%sigvz
         vdx = 0.0; vdy = 0.0; vdz = this%gamma
         sigx = this%sigx; sigy = this%sigy; sigz = this%sigz
         x0 = this%bcx; y0 = this%bcy; z0 = this%bcz
         cx = (/this%cx1,this%cx2,this%cx3/); cy = (/this%cy1,this%cy2,this%cy3/)
         lquiet = this%quiet
         idimp = size(part3d,1); npmax = size(part3d,2)
         noff = fd%getnoff()
         edges(1) = noff(1); edges(3) = noff(2)
         edges(2) = edges(1) + fd%getnd2p()
         edges(4) = edges(3) + fd%getnd3p()         
         
         call PRVDIST32_RANDOM(pt,this%qm,edges,npp,nps,vtx,vty,vtz,vdx,vdy,&
         &vdz,npx,npy,npz,nx,ny,nz,ipbc,idimp,npmax,1,1,4,sigx,sigy,sigz,&
         &x0,y0,z0,cx,cy,lquiet,ierr)

         if (ierr /= 0) then
            write (erstr,*) 'PRVDIST32_RANDOM error'
            call this%err%equit(class//sname//erstr)
         endif
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine dist3d_000
!
      subroutine init_fdist3d_001(this,input,i)
      
         implicit none
         
         class(fdist3d_001), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i        
! local data
         integer :: npf,npx,npy,npz,npmax
         real :: qm,sigx,sigy,bcx,bcy,bcz,sigvx,sigvy,sigvz
         real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np
         logical :: quiet
         real :: min, max, cwp, n0
         real :: alx, aly, alz, dx, dy, dz
         integer :: indx, indy, indz
         integer :: ii
         real :: sumz         
         character(len=20) :: sn,s1
         character(len=18), save :: sname = 'init_fdist3d_001:'
         
         this%sp => input%sp
         this%err => input%err
         this%p => input%pp

         call this%err%werrfl2(class//sname//' started')

         write (sn,'(I3.3)') i
         s1 = 'beam('//trim(sn)//')'

         call input%get('simulation.n0',n0)
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get('simulation.indz',indz)

         cwp=5.32150254*1e9/sqrt(n0)
         call input%get('simulation.box.x(1)',min)
         call input%get('simulation.box.x(2)',max)
         call input%get(trim(s1)//'.center(1)',bcx)
         bcx = bcx - min
         alx = (max-min) 
         dx=alx/real(2**indx)
         call input%get('simulation.box.y(1)',min)
         call input%get('simulation.box.y(2)',max)
         call input%get(trim(s1)//'.center(2)',bcy)
         bcy = bcy -min
         aly = (max-min) 
         dy=aly/real(2**indy)
         call input%get('simulation.box.z(1)',min)
         call input%get('simulation.box.z(2)',max)
         call input%get(trim(s1)//'.center(3)',bcz)
         bcz = bcz -min
         alz = (max-min) 
         dz=alz/real(2**indz)

         call input%get(trim(s1)//'.profile',npf)
         call input%get(trim(s1)//'.np(1)',npx)
         call input%get(trim(s1)//'.np(2)',npy)
         call input%get(trim(s1)//'.np(3)',npz)
         call input%get(trim(s1)//'.q',qm)
         call input%get(trim(s1)//'.sigma(1)',sigx)
         call input%get(trim(s1)//'.sigma(2)',sigy)
         call input%get(trim(s1)//'.sigma_v(1)',sigvx)
         call input%get(trim(s1)//'.sigma_v(2)',sigvy)
         call input%get(trim(s1)//'.sigma_v(3)',sigvz)
         call input%get(trim(s1)//'.centroid_x(1)',cx1)
         call input%get(trim(s1)//'.centroid_x(2)',cx2)
         call input%get(trim(s1)//'.centroid_x(3)',cx3)
         call input%get(trim(s1)//'.centroid_y(1)',cy1)
         call input%get(trim(s1)//'.centroid_y(2)',cy2)
         call input%get(trim(s1)//'.centroid_y(3)',cy3)
         call input%get(trim(s1)//'.quiet_start',quiet)
         call input%get(trim(s1)//'.gamma',gamma)
         call input%get(trim(s1)//'.peak_density',np)
         call input%get(trim(s1)//'.npmax',npmax)
         call input%get(trim(s1)//'.piecewise_fz',this%fz)
         call input%get(trim(s1)//'.piecewise_z',this%z)

         sumz = 0.0
         do ii = 2, size(this%z)
            if (this%z(ii)<=this%z(ii-1)) then
               write (erstr,*) 'Piecewise_z is not monotonically increasing'
               call this%err%equit(class//sname//erstr)
               return
            end if
            sumz = sumz + (this%fz(ii) + this%fz(ii-1))*&
            &(this%z(ii) - this%z(ii-1))*0.5
         end do
         this%z = this%z/dz
         this%npf = npf
         this%npx = npx
         this%npy = npy
         this%npz = npz
         this%npmax = npmax
         qm = qm/abs(qm)*abs(np)*(2*3.1415926535897932)*sigx*sigy*sumz
         qm = qm/dx/dy/dz
         qm = qm/npx
         qm = qm/npy
         qm = qm/npz
         this%qm = qm
         this%bcx = bcx/dx
         this%bcy = bcy/dy
         this%bcz = bcz/dz
         this%sigx = sigx/dx
         this%sigy = sigy/dy
         this%sigvx = sigvx
         this%sigvy = sigvy
         this%sigvz = sigvz
         this%cx1 = cx1*dz*dz/dx
         this%cx2 = cx2*dz/dx
         this%cx3 = cx3/dx
         this%cy1 = cy1*dz*dz/dy
         this%cy2 = cy2*dz/dy
         this%cy3 = cy3/dy
         this%gamma = gamma
         this%np = np
         this%quiet = quiet

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist3d_001
!
      subroutine dist3d_001(this,part3d,npp,fd)
      
         implicit none
         
         class(fdist3d_001), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part3d
         integer, intent(inout) :: npp
         class(ufield3d), intent(in), pointer :: fd
! local data
! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
         real, dimension(:,:), pointer :: pt => null()
         integer :: npx, npy, npz, nx, ny, nz, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real :: sigx, sigy, sigz, x0, y0, z0
         real, dimension(3) :: cx, cy
         real, dimension(4) :: edges
         integer, dimension(2) :: noff
         integer :: nps=1
         logical :: lquiet = .false.
         integer :: idimp, npmax, ierr = 0
         integer :: nzf, i, j
         real, dimension(:), allocatable :: zf
         character(len=18), save :: sname = 'dist3d_001:'

         call this%err%werrfl2(class//sname//' started')
         
         npx = this%npx; npy = this%npy; npz = this%npz
         nzf = size(this%z)
         nx = fd%getnd1(); ny = fd%getnd2(); nz = fd%getnd3()
         if ((this%z(1)>=nz) .or. (this%z(nzf)<=0)) then
            npp = 0
            return
         end if
         ipbc = this%sp%getpsolver()
         pt => part3d
         vtx = this%sigvx; vty = this%sigvy; vtz = this%sigvz
         vdx = 0.0; vdy = 0.0; vdz = this%gamma
         sigx = this%sigx; sigy = this%sigy
         x0 = this%bcx; y0 = this%bcy; z0 = this%bcz
         cx = (/this%cx1,this%cx2,this%cx3/); cy = (/this%cy1,this%cy2,this%cy3/)
         lquiet = this%quiet
         idimp = size(part3d,1); npmax = size(part3d,2)
         noff = fd%getnoff()
         edges(1) = noff(1); edges(3) = noff(2)
         edges(2) = edges(1) + fd%getnd2p()
         edges(4) = edges(3) + fd%getnd3p()
         allocate(zf(nz))
         zf = 0.0   
         do i = 1, nz
            do j = 2, nzf
               if ((i>=this%z(j-1)) .and. (i<this%z(j))) then
                  zf(i) = this%fz(j) + (this%fz(j-1)-this%fz(j))/&
                  &(this%z(j-1)-this%z(j))*(real(i)-this%z(j))
                  exit
               end if
            end do
         end do
         
         call PRVDIST32_RAN_PFL(pt,this%qm,edges,npp,nps,x0,y0,z0,sigx,sigy,&
         &vtx,vty,vtz,vdx,vdy,vdz,cx,cy,npx,npy,npz,nx,ny,nz,ipbc,idimp,&
         &npmax,1,1,4,zf,lquiet,ierr)

         if (ierr /= 0) then
            write (erstr,*) 'PRVDIST32_RAN_PFL error'
            call this%err%equit(class//sname//erstr)
         endif
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine dist3d_001
!
      subroutine init_fdist3d_002(this,input,i)
      
         implicit none
         
         class(fdist3d_002), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i        
! local data
         integer :: npf,npx,npy,npz,npmax
         real :: qm,sigx,sigy,sigz,bcx,bcy,bcz,sigvz
         real, dimension(2) :: alpha, beta, emit
         real :: gamma,np
         logical :: quiet, evol
         real :: min, max, cwp, n0
         real :: alx, aly, alz, dx, dy, dz
         integer :: indx, indy, indz         
         character(len=20) :: sn,s1
         character(len=18), save :: sname = 'init_fdist3d_002:'
         
         this%sp => input%sp
         this%err => input%err
         this%p => input%pp

         call this%err%werrfl2(class//sname//' started')

         write (sn,'(I3.3)') i
         s1 = 'beam('//trim(sn)//')'

         call input%get('simulation.n0',n0)
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get('simulation.indz',indz)

         cwp=5.32150254*1e9/sqrt(n0)
         call input%get('simulation.box.x(1)',min)
         call input%get('simulation.box.x(2)',max)
         call input%get(trim(s1)//'.center(1)',bcx)
         bcx = bcx - min
         alx = (max-min) 
         dx=alx/real(2**indx)
         call input%get('simulation.box.y(1)',min)
         call input%get('simulation.box.y(2)',max)
         call input%get(trim(s1)//'.center(2)',bcy)
         bcy = bcy -min
         aly = (max-min) 
         dy=aly/real(2**indy)
         call input%get('simulation.box.z(1)',min)
         call input%get('simulation.box.z(2)',max)
         call input%get(trim(s1)//'.center(3)',bcz)
         bcz = bcz -min
         alz = (max-min) 
         dz=alz/real(2**indz)

         call input%get(trim(s1)//'.profile',npf)
         call input%get(trim(s1)//'.np(1)',npx)
         call input%get(trim(s1)//'.np(2)',npy)
         call input%get(trim(s1)//'.np(3)',npz)
         call input%get(trim(s1)//'.q',qm)
         call input%get(trim(s1)//'.sigmaz',sigz)
         call input%get(trim(s1)//'.sigma_vz',sigvz)
         call input%get(trim(s1)//'.alpha(1)',alpha(1))
         call input%get(trim(s1)//'.alpha(2)',alpha(2))
         call input%get(trim(s1)//'.beta(1)',beta(1))
         call input%get(trim(s1)//'.beta(2)',beta(2))
         call input%get(trim(s1)//'.emittance(1)',emit(1))
         call input%get(trim(s1)//'.emittance(2)',emit(2))
         call input%get(trim(s1)//'.quiet_start',quiet)
         call input%get(trim(s1)//'.gamma',gamma)
         call input%get(trim(s1)//'.peak_density',np)
         call input%get(trim(s1)//'.npmax',npmax)
         call input%get(trim(s1)//'.evolution',evol)

         this%npf = npf
         this%npx = npx
         this%npy = npy
         this%npz = npz
         this%npmax = npmax
         sigx = sqrt(beta(1)*emit(1)/gamma)
         sigy = sqrt(beta(2)*emit(2)/gamma)
         qm = qm/abs(qm)*abs(np)*(2*3.1415926535897932)**1.5*sigx*sigy*sigz
         qm = qm/dx/dy/dz
         qm = qm/npx
         qm = qm/npy
         qm = qm/npz
         this%qm = qm
         this%bcx = bcx/dx
         this%bcy = bcy/dy
         this%bcz = bcz/dz
         this%sigz = sigz/dz
         this%sigvz = sigvz
         this%alphax = alpha(1)
         this%alphay = alpha(2)
         this%betax = beta(1)/dx
         this%betay = beta(2)/dy
         this%emitx = emit(1)/dx
         this%emity = emit(2)/dy
         this%gamma = gamma
         this%np = np
         this%quiet = quiet
         this%evol = evol

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist3d_002
!
      subroutine dist3d_002(this,part3d,npp,fd)
      
         implicit none
         
         class(fdist3d_002), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part3d
         integer, intent(inout) :: npp
         class(ufield3d), intent(in), pointer :: fd
! local data1

! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
         real, dimension(:,:), pointer :: pt => null()
         integer :: npx, npy, npz, nx, ny, nz, ipbc
         real :: vtz, vdx, vdy, vdz
         real :: sigz, x0, y0, z0
         real :: alphax, betax, emitx
         real :: alphay, betay, emity
         real, dimension(4) :: edges
         integer, dimension(2) :: noff
         integer :: nps=1
         logical :: lquiet = .false.
         integer :: idimp, npmax, ierr = 0
         character(len=18), save :: sname = 'dist3d_002:'

         call this%err%werrfl2(class//sname//' started')
         
         npx = this%npx; npy = this%npy; npz = this%npz
         nx = fd%getnd1(); ny = fd%getnd2(); nz = fd%getnd3()
         ipbc = this%sp%getpsolver()
         pt => part3d
         sigz = this%sigz; vtz = this%sigvz
         vdx = 0.0; vdy = 0.0; vdz = this%gamma
         alphax = this%alphax; alphay = this%alphay
         betax = this%betax; betay = this%betay
         emitx = this%emitx; emity = this%emity
         x0 = this%bcx; y0 = this%bcy; z0 = this%bcz
         lquiet = this%quiet
         idimp = size(part3d,1); npmax = size(part3d,2)
         noff = fd%getnoff()
         edges(1) = noff(1); edges(3) = noff(2)
         edges(2) = edges(1) + fd%getnd2p()
         edges(4) = edges(3) + fd%getnd3p()         
         
         call PRVDIST32_TWISS(pt,this%qm,edges,npp, nps, alphax,alphay,   &
         &betax,betay, emitx,emity,sigz,vdx,vdy,vdz,vtz,npx,npy,npz,idimp,&
         &npmax,nx,ny,nz,x0,y0,z0,1,1,4,ierr,vdz,lquiet)

         ! call PRVDIST32_RANDOM(pt,this%qm,edges,npp,nps,vtx,vty,vtz,vdx,vdy,&
         ! &vdz,npx,npy,npz,nx,ny,nz,ipbc,idimp,npmax,1,1,4,sigx,sigy,sigz,&
         ! &x0,y0,z0,cx,cy,lquiet,ierr)

         if (ierr /= 0) then
            write (erstr,*) 'PRVDIST32_RANDOM error'
            call this%err%equit(class//sname//erstr)
         endif
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine dist3d_002
!
      subroutine init_fdist3d_003(this,input,i)

         implicit none

         class(fdist3d_003), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i
! local data
         integer :: npt, npmax, npf
         real :: bcx, bcy, bcz, tc, qm
         logical :: evol
         real :: min, max, cwp, n0
         real :: alx, aly, alz, dx, dy, dz
         integer :: indx, indy, indz
         character(len=20) :: sn,s1
         character(len=18), save :: sname = 'init_fdist3d_003:'

         this%sp => input%sp
         this%err => input%err
         this%p => input%pp

         call this%err%werrfl2(class//sname//' started')

         write (sn,'(I3.3)') i
         s1 = 'beam('//trim(sn)//')'

         call input%get('simulation.n0',n0)
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get('simulation.indz',indz)

         cwp=5.32150254*1e9/sqrt(n0)
         this%cwp = cwp
         call input%get('simulation.box.x(1)',min)
         call input%get('simulation.box.x(2)',max)
         call input%get(trim(s1)//'.center(1)',bcx) !beam center
         bcx = bcx - min
         alx = (max-min) 
         dx=alx/real(2**indx)
         call input%get('simulation.box.y(1)',min)
         call input%get('simulation.box.y(2)',max)
         call input%get(trim(s1)//'.center(2)',bcy)
         bcy = bcy -min
         aly = (max-min) 
         dy=aly/real(2**indy)
         call input%get('simulation.box.z(1)',min)
         call input%get('simulation.box.z(2)',max)
         call input%get(trim(s1)//'.center(3)',bcz)
         bcz = bcz -min
         alz = (max-min) 
         dz=alz/real(2**indz)

         call input%get(trim(s1)//'.profile',npf)
         call input%get(trim(s1)//'.np',npt) !npt多少粒子
         call input%get(trim(s1)//'.npmax',npmax)!与并行核有关，动态分配内存会产生内存碎片，程序运算很慢/一开始内存取好，尽量去大一点
         call input%get(trim(s1)//'.evolution',evol)!beam shape等不变，粒子不变
         call input%get(trim(s1)//'.file_name',this%file)
         call input%get(trim(s1)//'.total_particle_num',tc)
         call input%get(trim(s1)//'.q',qm)


         this%npf = npf
         this%npt = npt
         this%npmax = npmax
         this%bcx = bcx/dx !cell size归一化
         this%bcy = bcy/dy
         this%bcz = bcz/dz
         this%dx = dx
         this%dy = dy
         this%dz = dz
         this%evol = evol
         this%qm = qm*abs(tc)*(1.0e12)/(cwp**3)/n0/dx/dy/dz/npt

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist3d_003
!
      subroutine dist3d_003(this,part3d,npp,fd)!子程序，可改变调用参数值、全局变量的值，可有任意返回值
      
         implicit none ! Require all variables to be explicitly declared
         
         class(fdist3d_003), intent(inout) :: this !class关键词实现动态绑定，inout提供/返回数据
         real, dimension(:,:), pointer, intent(inout) :: part3d !part3d可指向实型目标变量
         integer, intent(inout) :: npp
         class(ufield3d), intent(in), pointer :: fd !输入参数
! local data1

! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
         real, dimension(:,:), pointer :: part
         integer :: npt, nx, ny, nz, ipbc, i !npt? ipbc?
         real :: x0, y0, z0, dx, dy, dz, cwp, qm !c/wp;xyz[um];vxvyvz[gamma*v/c]
         real, dimension(4) :: edges
         integer, dimension(2) :: noff !?
         integer :: nps !?
         integer :: idimp, npmax, ierr !idimp:
         real tempx, tempy, tempz, tempvx, tempvy, tempvz, tempq !电荷密度
         character(len=18), save :: sname = 'dist3d_003:' !len最大长度

         call this%err%werrfl2(class//sname//' started') !werrfl2? 合并

         nps = 1 
         ierr = 0
         npt = this%npt !total number
         nx = fd%getnd1(); ny = fd%getnd2(); nz = fd%getnd3()
         ipbc = this%sp%getpsolver() !?
         part => part3d 
         x0 = this%bcx; y0 = this%bcy; z0 = this%bcz !beam center/box center?
         dx = this%dx; dy = this%dy; dz = this%dz; cwp = this%cwp
         qm = this%qm
         idimp = size(part3d,1); npmax = size(part3d,2)
         noff = fd%getnoff()
         edges(1) = noff(1); edges(3) = noff(2)
         edges(2) = edges(1) + fd%getnd2p()
         edges(4) = edges(3) + fd%getnd3p()         
         
         open(unit=85,file=trim(this%file),status='old',action='read') !返回字符串尾端多余空格清除后的字符串
         
         do i = 1, npt
         
            read(85,*,iostat=ierr) tempz,tempx,tempy,tempvz,tempvx,tempvy,tempq !isotat的值用于判断输入成功否。
            
            tempx = tempx/cwp/dx + x0
            tempy = tempy/cwp/dy + y0
            tempz = -tempz/cwp/dz + z0 !z=ct-z

            if ((tempx>=1) .and. tempx<=(nx-1) .and. (tempy >= edges(1))& !是用来判断粒子是否属于当前这个分区的。
            & .and. (tempy < edges(2)) .and. (tempz >= edges(3)) .and. (t& !因为是并行程序，每个核负责空间里面的一块儿
            &empz < edges(4))) then !读进来一个粒子，要判断这个粒子是否位于当前这个核所负责的区域，在就留下，不在就读下一个粒子。
              part(1,nps) = tempx
              part(4,nps) = tempvx
              part(2,nps) = tempy
              part(5,nps) = tempvy
              part(3,nps) = tempz
              part(6,nps) = tempvz
              part(7,nps) = tempq*qm
              nps = nps + 1
              cycle
            endif     
         enddo
         
         npp = nps - 1 !当前核粒子数
         close(85)

         if (ierr /= 0) then
            write (erstr,*) 'Read beam file error'
            call this%err%equit(class//sname//erstr)
         endif
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine dist3d_003
!
      subroutine init_fdist3d_100(this,input,i)
      
         implicit none
         
         class(fdist3d_100), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i        
! local data
         integer :: xppc,yppc,zppc,npmax,npf
         real :: qm, bcx, bcy, bcz
         real, dimension(3) :: lb, ub
         real :: r1, r2
         real :: gamma, np
         real :: min, max, cwp, n0
         real :: alx, aly, alz, dx, dy, dz
         integer :: indx, indy, indz         
         character(len=20) :: sn,s1
         character(len=18), save :: sname = 'init_fdist3d_100:'
         
         this%sp => input%sp
         this%err => input%err
         this%p => input%pp

         call this%err%werrfl2(class//sname//' started')

         write (sn,'(I3.3)') i
         s1 = 'beam('//trim(sn)//')'

         call input%get('simulation.n0',n0)
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get('simulation.indz',indz)

         cwp=5.32150254*1e9/sqrt(n0)
         call input%get('simulation.box.x(1)',min)
         call input%get('simulation.box.x(2)',max)
         call input%get(trim(s1)//'.center(1)',bcx)
         call input%get(trim(s1)//'.range_x(1)',lb(1))
         call input%get(trim(s1)//'.range_x(2)',ub(1))
         lb(1) = lb(1) - min
         ub(1) = ub(1) - min
         bcx = bcx - min
         alx = (max-min) 
         dx=alx/real(2**indx)
         call input%get('simulation.box.y(1)',min)
         call input%get('simulation.box.y(2)',max)
         call input%get(trim(s1)//'.center(2)',bcy)
         call input%get(trim(s1)//'.range_y(1)',lb(2))
         call input%get(trim(s1)//'.range_y(2)',ub(2))
         lb(2) = lb(2) - min
         ub(2) = ub(2) - min
         bcy = bcy -min
         aly = (max-min) 
         dy=aly/real(2**indy)
         call input%get('simulation.box.z(1)',min)
         call input%get('simulation.box.z(2)',max)
         call input%get(trim(s1)//'.center(3)',bcz)
         call input%get(trim(s1)//'.range_z(1)',lb(3))
         call input%get(trim(s1)//'.range_z(2)',ub(3))
         lb(3) = lb(3) - min
         ub(3) = ub(3) - min
         bcz = bcz -min
         alz = (max-min) 
         dz=alz/real(2**indz)

         call input%get(trim(s1)//'.profile',npf)
         call input%get(trim(s1)//'.ppc(1)',xppc)
         call input%get(trim(s1)//'.ppc(2)',yppc)
         call input%get(trim(s1)//'.ppc(3)',zppc)
         call input%get(trim(s1)//'.q',qm)
         call input%get(trim(s1)//'.np',np)
         call input%get(trim(s1)//'.gamma',gamma)
         call input%get(trim(s1)//'.npmax',npmax)
         call input%get(trim(s1)//'.r1',r1)
         call input%get(trim(s1)//'.r2',r2)

         this%npf = npf
         this%npmax = npmax
         this%xppc = xppc
         this%yppc = yppc
         this%zppc = zppc
         this%bcx = bcx/dx
         this%bcy = bcy/dy
         this%bcz = bcz/dz
         this%lb(1) = lb(1)/dx
         this%lb(2) = lb(2)/dy
         this%lb(3) = lb(3)/dz
         this%ub(1) = ub(1)/dx
         this%ub(2) = ub(2)/dy
         this%ub(3) = ub(3)/dz
         this%r1 = r1/dx
         this%r2 = r2/dx
         this%gamma = gamma
         this%np = np

         qm = qm/abs(qm)/xppc/yppc/zppc
         this%qm = qm

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist3d_100
!
      subroutine dist3d_100(this,part3d,npp,fd)
      
         implicit none
         
         class(fdist3d_100), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part3d
         integer, intent(inout) :: npp
         class(ufield3d), intent(in), pointer :: fd
! local data1
         real, dimension(:,:), pointer :: pt => null()
         integer :: xppc,yppc,zppc
         real :: qm, bcx, bcy, bcz
         real, dimension(3) :: lb, ub
         real :: r1, r2, tx, ty, tz, a1, a2, a3, a4
         real :: gamma, np
         integer, dimension(2) :: noff
         integer :: nps
         integer :: i, j, k, ix, iy, iz
         integer :: idimp, npmax, ierr
         character(len=18), save :: sname = 'dist3d_100:'

         call this%err%werrfl2(class//sname//' started')
         
         nps = 1
         ierr = 0
         xppc = this%xppc; yppc = this%yppc; zppc = this%zppc;
         noff = fd%getnoff()
         lb(1) = max(this%lb(1),0.0)
         ub(1) = min(this%ub(1),real(fd%getnd1()))
         if (this%lb(2) > real(noff(1)+fd%getnd2p()) .or. &
         &this%ub(2) < real(noff(1))) then
            npp = 0
            return
         else
            lb(2) = max(this%lb(2),real(noff(1)))
            ub(2) = min(this%ub(2),real(noff(1)+fd%getnd2p()))
         end if
         if (this%lb(3) > real(noff(2)+fd%getnd3p()) .or. &
         &this%ub(3) < real(noff(2))) then
            npp = 0
            return
         else
            lb(3) = max(this%lb(3),real(noff(2)))
            ub(3) = min(this%ub(3),real(noff(2)+fd%getnd3p()))
         end if

         pt => part3d
         r1 = this%r1; r2 = this%r2; np = this%np
         gamma = this%gamma
         bcx = this%bcx; bcy = this%bcy; bcz = this%bcz
         qm = this%qm
         bcx = this%bcx; bcy = this%bcy; bcz = this%bcz

         do i=int(lb(1)), int(ub(1))-1
            do j=int(lb(2)), int(ub(2))-1
               do k = int(lb(3)), int(ub(3))-1
                  do ix = 0, xppc-1
                     do iy= 0, yppc-1
                        do iz= 0, zppc-1
                           tx = (ix + 0.5)/xppc + i - 1
                           ty = (iy + 0.5)/yppc + j - 1
                           tz = (iz + 0.5)/zppc + k - 1
                           pt(1,nps) = tx
                           pt(2,nps) = ty
                           pt(3,nps) = tz
                           pt(4,nps) = 0.0
                           pt(5,nps) = 0.0
                           pt(6,nps) = gamma
                           a1 = abs(tz - bcz)
                           a2 = sqrt((tx-bcx)**2+(ty-bcy)**2)
                           a3 = abs(a2-r1)
                           a4 = a1**2+a3**2
                           pt(7,nps) = qm*np*exp(-a4/r2**2*0.5)
                           nps = nps + 1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo

         npp = nps - 1
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine dist3d_100
!
      end module fdist3d_class