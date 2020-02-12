c change from qro_v6a3.for   10/22/2018
c    archp(ilk,0:2,0:MXPS,icy,icm) is modified to archp(ilk,0:2,icy,icm),  
c    aflow(ili,0:2,0:MXPS,icy,icm) is modified to aflow(ilk,0:2,icy,icm),  
c    to produce output, ignored point source components, and only total flow is considerted
c    no modeling difference
c
c
C February, 2018, modified from \Chino11\FOR|ro_v6a3.FOR
C
C for water transfer example, see F:\Projects\CBRMPU16\Truing-Up\Opt23a\
C
C 10/31/2016 modified from ro_v6a2.for to handle water transfer between storm water basins
C
C Procedure to evaluation of pumping stormwater to upstream basins
C    a. run this or previous version v6a1, 
C       if nqls > 0, two output file name must be provided:
C       1) monthly percolation *.lmp and 2) monthly evaporation *.lme
C       if a thrid file is specified *.sto, end-of-day storage will be listed for all basins in nqls list
C       include all basins that water is diverted from and basins which will receive diversion
C    b. compare end of day storage with basin capacity, and storage at the diversion source basin, 
C       create daily pumping schedule file
C       Include daily delivery data to either boundary inflow or point discharge file with proper name
C    c. ccRerun with this version
C       modify node file to receive daily diversion
C       modify link file to specify diversion basin,  the ID number of predefined diversion for the basin on line for number data entry
C       modify con file to sprcify number of predefined diversion on the line for number of point discharge
C       if a foutth file name is specified for nqls basins, echo actual diversion will be printed for check.
C  Procedure a through c could be repeated with updated data to maximize.
C
C 7/21/2014 modified to print daily storage in recharge basins
c
c 11/2008 modified from ro_v4r4.for to include multipurpose basin operation
c
c 3/18/08 modified to facilitate PEST application
c
c 11/11/07 modification to include water quality
c
c 7/26/07 variable perc rate is allowed for recharge basins
c
c route_q6.for
c - this version each tracks point-sources in flow and constituent movement
c - reservoir routine is revised
c - known diversion can be specified as point-source with negative flow
C
c route_q5.for
c this version is modified from route_q4 to
c   summarize stream percolation for reaches specified in link file
c reservoir routine is revised
c
c route_q4.for
c this version is modified from route_n9 to include water quality simulation
c rxn(ilk)  - positive number, 0 < rxn <= 1, to be multiplied to incoming
c water quality for inflows can be specified in three methods:
c   1 daily value
c   2 flow-conc relationship
c   3 constant value
c the method can be spefied for the source of water, runoff model output,
c boundary inflow, or point source
c the method is spefied at the top of the file
c
c Modified by Jeff on 10/01
c         This version includes options for:
c            Node, link, HA names are 8 alphanumeric chracters
c            separate file for recorded boundary inflow file
c            separate file for wastewater discharge file
c  

C Files opened
c  Unit  File description
!    1    main input file
!    2    link file
!    3    node file
!    4    open channel parameters
!         open channel flow-wetted surface area rating table
!    5    daily evaporation data file
!    6    main output file
!    7    storm prediction rainfall data file
!    8
!    9
!   10    diversion output file (optional)   used for debugging purpose
!   11    runoff model output flow data file
!   12    boundary inflow data file
!   13    point source flow data file 
!   14    HSA character file, share file with RF module,  CNzone, Rain zone, GWM fraction, and LID onsite removal
!   15    predefined diversion file, if npdd>0  (daily diversion data)
!   16    diversion scheme spec file, if file name is given after *.STO or *.ADV
!   21    Runoff model output, TDS on daily runoff from each HSA
!   22    boundary inflow TDS data file, either daily TDS, fixed value, or rating table
!   23    point source TDS data file 
!   31    Runoff model output, TIN on daily runoff from each HSA
!   32    boundary inflow TIN data file, either daily TIN, fixed value, or rating table
!   33    point source TIN data file
!   61    output file to print flow at user specified link or node, daily, nplk
!   62    output file to print TDS  at user specified link or node
!   63    output file to print TIN  at user specified link or node
!   64    output file to print flow at user spefied node (for storm water and each point source), monthly, npnd
!   65    output file to print TDS  at user spefied node (for storm water and each point source)
!   66    output file to print TIN  at user spefied node (for storm water and each point source)
!   67    output file to print percolation in user spedifed reaches (for storm water and each point source), month, nprch
!   68    output file to print TDS in percolation in user spedifed reaches (for storm water and each point source)
!   69    output file to print TIN in percolation in user spedifed reaches (for storm water and each point source)
!   81    *.LMP --- monthly percolation output for recharge basins, nqls
!   82    *.LMT --- TDS in monthly percolation output for recharge basins
!   83    *.LMN --- NO3 in monthly percolation output for recharge basins
!   84    *.LME --- monthly evaporation output for recharge basins
!   85    *.STO --- end of day storage for nqls basins
!   86    *.STT --- TDS
!   86    *.STN --- TIN

c   91    *.ADV --- actual diversion for npdd reservoirs, if npdd>0    !v6a2
c   92    *.DSH --- actual daily pumping for specified diversion scheme

c  TDS and TIN data input option
c     mopt(i,j)
c     i = 1     precipitation runoff, nha   j = 1   TDS
c         2     boundary inflow, nbf            2   TIN
c         3     point source, nps
c     mopt = 1  daily TDS and TIN concentration is provided, --- Runoff module should generate TDS and TIN files
c            2  constant concentration, for i=2, for planning purpose, used for WLAM scenarios
c            3  flow (inches) vs concentration table, for i=1,2
c            0 or blank is considered as 1

      include 'RO_DIM_3.MAX'

      logical      okfile

      character*50 file0,lfile,nfile,hafile,bifile,psfile,
     &             efile,ofile,divfile,bl50,chfile,
     &             ffile,lidfile,czfile,pdfile,ifile,opfile(3,2),
     &             prfile        ! precipitation file as storm forecaster   !r5
      character*40 fmtin,rfmte,rfmtp
      integer      t
ccccc integer      t,mopt(3,2)
      character*8  ch8(10), bl8
      real         rtemp(MXHA),rtemp2(MXHA)     ! temperary variables

      real         alpha(10),dperc(10), dmin(10),
     &             pcap(10), pmcap(10), prcap(10),
     &             ttrans(10), ptrans(10), pimp(10)
      integer      ipipe(10),lkatts(10), kratts(10),
     &                       lkatte(10), kratte(10)
      
      integer ntconc(14,2)
	  real    tconc(14,2,10,2)
	  
c---- general model control properties
      include 'GenParm.VAR'

c---- link/node properties
      include 'LKprpt2.var'

c---- ha property
      include 'HAprpt.var'

c---- reservoir property
      include 'RVprpt.VAR'

c---- diversion properties
      include 'DVprpt.VAR'

c---- table properties
      include 'TBprpt.VAR'

c---- output variables
      include 'OVprpt.VAR'

      dimension    evapor(MXEV)

      integer      kpql(MXRV)
      
      ! unit convertion from (cfs * day) to (acre-ft)
      cf1 = 86400./43560.
      ! unit conversion from (cfs * mg/l * day) to (tons)
      cf2 = 86400*28.32/(2000.*453600.)
      delimiter = '|'
      cstar = '*'

      ierror = .FALSE.
      bl8 = '        '  
      bl50 = '                                                  '
      velocity = 1.0

c     initialize variables
      call vrinit()

      archp = 0.
      aflow = 0.
      kpql = 0     !v6a1

c
c     read in control parameters from control file
c

      call getcon(file0,n1)

      open(1,file=file0,status='OLD')
      call skipline(1,'*')
      read(1,1003) newts,nhat,nbf
	  read(1,'(2i8)') nps, npdd     ! number of point source, number pre-defined diversion
	  read(1,1003) idebug,iqopt
      nhat1 = nhat + nbf
      nhat2 = nhat1 + nps
      qopt = .FALSE.
      if(iqopt.ne.0.) qopt = .TRUE.
      ldebug = .FALSE.
      if(idebug.gt.0) ldebug = .TRUE.
      read(1,'(f8.0,i8)') rainff,nffday
      read(1,'(f8.0,i8)') rainss
      write(*,'(a,2f8.2,i8)') ' Rainss, rainff, nffday =',
     &                        rainss,rainff,nffday
         ! rainff, minimum rainfall for first flush
         ! ndffday, number of dry days to qualify as first flush
         ! rainss, minimum rainfall for significant storm

 1003 format(i8)


      
      call OpnIF2(1, 2, lfile,cstar)         ! link definition file
      call OpnIF2(1, 3, nfile,cstar)         ! node definition file
      call OpnIF2(1, 4, ffile,cstar)         ! reach parameter file

      read(4,1003) ncnz     ! data used by Runoff module
      do i=1,ncnz
        read(4,1003) j
      end do
      call skipline(4,cstar)
      read(4,1003) nrainz   ! data used by Runoff module
      do i=1,nrainz
        read(4,1003) j
      end do
      call skipline(4,cstar)
      read(4,1003) nrchz    !
      do i=1,nrchz
        read(4,'(i8,4f8.0)') j,(rvrch(j,k),k=1,3)
        if(i.ne.j) then
          write(*,'(a)') ' Error in parameter file'
          stop
        end if
      end do
      close(4)


! flow input files
      call OpnIF2(1,11,hafile,cstar)         ! runoff model output file
      call OpnIF2(1,12,bifile,cstar)         ! recorded boundary inflow file
      call OpnIF2(1,13,psfile,cstar)         ! recorded point source file
      
!****************************************************************
c
c     read hydrologic area names from flow data files
c

      read(11,'(8x,25a8)') (haname(kk),kk=1,nhat)
      read(11,'((8x,25f8.0))') (haarea(kk),kk=1,nhat)

      ! read start date
      read(11,'(i4,2i2)') isyr,ismo,isdy
      backspace(11)
      write(*,'(a,i4.4,2i2.2)') ' Simulation starts at ',isyr,ismo,isdy
      
      if(nbf.gt.0) then
        read(12,'(8x,30a8)')
     &            (haname(kk),kk=nhat+1,nhat1)
        ! position the file at the start of simulation
        CALL LocTS1(12,isyr,ismo,isdy)
      end if
      continue

      if(nps.gt.0) then
        read(13,'(8x,50a8)')
     &            (haname(kk),kk=nhat1+1,nhat2)
        ! position the file at the start of simulation
        CALL LocTS1(13,isyr,ismo,isdy)
      end if
      continue

      call OpnIF2(1,14,czfile,cstar)         ! ha area parameters
      read(14,'(a40)') fmtin
      do j=1,nhat
        read(14,fmtin,err=10,end=10) ii,ch8(1),                       
     &             areaha, i1,i2,r1,r2      ! CN zone, rain zone, GWM area fraction(r1) used by Runoff module
        ! find HAID
        flid(j) = r2                     ! flow removed on site by LID
        !this option is not available yet
      end do
 10   close(14)

	  if(npdd.gt.0) call OpnIF2(1,15,pdfile,cstar)

      if(qopt) then
        ! TDS input files
        call OpnIF2(1,21,opfile(1,1),cstar)     ! RF model runoff TDS
        if(opfile(1,1).ne.bl50) read(21,1003) mopt(1,1)
        if(mopt(1,1).eq.0) mopt(1,1) = 1
        ! only option 1 is available for storm runoff water quality
        if(mopt(1,1).gt.1) stop
        call LocTS1(21,isyr,ismo,isdy)
 
        call OpnIF2(1,22,opfile(2,1),cstar)     ! Boundary inflow TDS
        if(opfile(2,1).ne.bl50) read(22,1003) mopt(2,1)
        if(mopt(2,1).eq.0) mopt(2,1) = 1
        if(mopt(2,1).le.1) then                 ! daily value
          call LocTS1(22,isyr,ismo,isdy)
        elseif(mopt(2,1).eq.2) then             ! constant TDS
          read(22,'(8x,50f8.0)') (cfixed(i,1),i=nhat+1,nhat1)
        else                                    ! flow vs TDS
          read(22,'(8x,30i8)') (ptrc(i,1),i=nhat+1,nhat1)  ! table number
          read(22,'(8x,30f8.0)') (flowconv(i,1),i=nhat+1,nhat1)
ccc          do while (.true.)
ccc            read(22,'(2i8)',err=20, end=20) ibf,ibfen
ccc            if(ibf.eq.0) goto 20
ccc            ntabq(ibf,1) = ibfen  ! number of entry for flow - tds table i
ccc            do i=1,ibfen
ccc              read(22,'(2f8.0)') fqtab(ii,1,i,1),fqtab(ii,1,i,2)
ccc                                 ! flow in inches, TDS
ccc            end do
ccc          end do            
ccc 20       continue
        end if  !mopt(2,1)
 
        call OpnIF2(1,23,opfile(3,1),cstar)     ! Point discharge TDS
        if(opfile(3,1).ne.bl50) read(23,1003) mopt(3,1)
        if(mopt(3,1).eq.0) mopt(3,1) = 1
        ! only option types 1 and 2 are allowed
        if(mopt(3,1).eq.3) stop
        if(mopt(3,1).le.1) call LocTS1(23,isyr,ismo,isdy)
        ! if option 2, read constant concentration
        if(mopt(3,1).eq.2) read(23,'(8x,50f8.0)')
     &                              (cfixed(i,1),i=nhat1+1,nhat2)
 
        ! TIN input files
        call OpnIF2(1,31,opfile(1,2),cstar)     ! RF model runoff TIN
        if(opfile(1,2).ne.bl50) read(31,1003) mopt(1,2)
        if(mopt(1,2).eq.0) mopt(1,2) = 1
        ! only option 1 is available
        if(mopt(1,1).gt.1) stop
        CALL LocTS1(31,isyr,ismo,isdy)
 
        call OpnIF2(1,32,opfile(2,2),cstar)     ! Boundary inflow TIN
        if(opfile(2,2).ne.bl50) read(32,1003) mopt(2,2)
        if(mopt(2,2).eq.0) mopt(2,2) = 1
        if(mopt(2,2).eq.1) then                 ! daily value
          call LocTS1(32,isyr,ismo,isdy)
        elseif(mopt(2,2).eq.2) then             ! constant TIN
          read(32,'(8x,50f8.0)') (cfixed(i,2),i=nhat+1,nhat1)
        else                                    ! flow vs TIN
          read(32,'(8x,30i8)') (ptrc(i,2),i=nhat+1,nhat1)
          read(32,'(8x,30f8.0)') (flowconv(i,2),i=nhat+1,nhat1)
ccc          do while (.true.)
ccc            read(22,'(2i8)',err=30, end=30) ibf,ibfen
ccc            if(ibf.eq.0) goto 30
ccc            ntabq(ibf,3) = ibfen  ! number of entry for flow - tds table i
ccc            do i=1,ibfen
ccc              read(22,'(2f8.0)') fqtab(ii,2,i,1),fqtab(ii,2,i,2)
ccc                                 ! flow in inches, TIN
ccc            end do
ccc          end do            
ccc 30       continue
        end if  ! mopt(2,2)

        call OpnIF2(1,33,opfile(3,2),cstar)     ! Point discharge TIN
        if(opfile(3,2).ne.bl50) read(33,1003) mopt(3,2)
        if(mopt(3,2).eq.0) mopt(3,2) = 1
        ! only option types 1 and 2 are allowed
        if(mopt(3,2).eq.3) stop
        if(mopt(3,2).le.1) CALL LocTS1(33,isyr,ismo,isdy)
        ! if option 2, read constant concentration
        if(mopt(3,2).eq.2) read(33,'(8x,50f8.0)')
     &                              (cfixed(i,2),i=nhat1+1,nhat2)
      end if  ! Qopt
      

! channel rating curves

      call OpnIF2(1, 4,chfile,cstar)         ! channel flow-width rating file
      nch = 0
      if(chfile.ne.bl50) then
        read(4,'(i8)') nch
        if(nch.gt.MXCH) then
          write(*,'(a)') ' Increase MXCH'
		  write(*,*) nch,MXCH
          stop
        end if
        do ich=1,nch
          read(4,'(a8,i6,50f6.0)') chname(ich),ntabf(ich),
     &                             (fwtab(ich,k,1),k=1,50)   ! flow
          read(4,'(8x,6x,50f6.0)') (fwtab(ich,k,2),k=1,50)   ! width
          call ctrim8r(chname(ich))
        end do
        close(4)
      end if

      if(qopt) then
        do ich=1,2
          call OpnIF2(1,9,ifile,cstar)                ! runoff-TDS or TIN tables
          do ilu=1,14   !!!!!!!!!!!!!!!!!!!!nlu
            read (9,'(2i8)') i,j
            if(i.ne.ilu) stop
            ntconc(ilu,ich) = j
            do i=1,j
              read(9,'(2f8.0)') tconc(ilu,ich,i,1),tconc(ilu,ich,i,2)
            end do
          end do
          close(9)
        end do
	  end if
	  
      call OpnIF2(1,5,efile,cstar)
      call LocTS1(5,isyr,ismo,isdy)
      call OpnIF2(1,7,prfile,cstar)
      call LocTS1(7,isyr,ismo,isdy)
      
      read(1,'(a50,i2)') ofile
        open(6,file=ofile,status='UNKNOWN')
        write(6,'(a,a)') ' Main Input File = ',file0

      read(1,'(a50)') divfile
      write(*,'(1x,a)') divfile
      if(divfile.ne.bl50) open(10,file=divfile,status='UNKNOWN')

      read(1,'(a50)') lidfile
      write(*,'(1x,a)') lidfile
      if(lidfile.ne.bl50) then
        open(14,file=lidfile,status='UNKNOWN')
        write(14,'(a8      ,25a8  ,/,(8x,25a8  ))') 'LID Flow',
     &                                  (haname(iha),iha=1,nhat)
      end if

      write(6,1003) newts,nhat,nbf,nps,idebug,iqopt
      ! if idebug = 0, no debug
      !             1, minimum debug
      !             2, full debug

 
      if(npdd.gt.0) then
        ! position the file at the start of new simulation
        call LocTS1(15,isyr,ismo,isdy) 
      end if

      do kk=1,nhat2
        call ctrim8r(haname(kk))
        write(6,'(i5,2x,a8,f8.0)') kk,haname(kk), haarea(kk)
        ! check for haname duplication
        if(kk.gt.1) then
          do jj=1,kk-1
            if(haname(jj).eq.haname(kk)) then
              write(*,'(a,2(i5,1x,a8))') ' Duplicate HANAME',
     &                                 jj,haname(jj),kk,haname(kk)
              ierror = .true.
            end if
          end do
        end if
      end do
      write(6,'(a1)') ' '

      call LNSetup()

      call rdoutp(file0)
      
      do i=1,nqls
        do k=1,nres
          if(lkatrv(k).eq.ipql(i)) then
            kpql(i) = k
            goto 50
          end if
        end do
        write(6,'(a,3i5)')' problem in nqls spec',i
		kpql(i) = 0    ! 1/29/19
		ipql(i) = 0    ! 1/29/19
 50     continue
      end do
      write(6,'(//,a)') ' Checking pointers'
      do i=1,nqls
         write(6,'(2i6,2x,a8,i6)')
     &                     i,ipql(i),lkname(ipql(i)),kpql(i)
      end do


      write(6,'(//,a)') 'Channel Flow-Width relationship'
      do ich=0,nch
        mtab = ntabf(ich)
        if(mtab.eq.0) goto 60
        write(6,'(i10,3(i10,10x))') ich, mtab
        do i=1,mtab
          write(6,'(i8,6f10.2)') i,(fwtab(ich,i,j),j=1,2)
        end do
 60     continue
      end do

      if(ierror) then
        write(*,'(a)') ' IERROR = .TRUE.'
        stop
      end if

      if(ldebug) then
        write(6,'(//,a)') 'Reservoir summary data'
        do k=1,nres
          write(6,'(a)') rvname(k)
          write(6,'(i5,a)') krvtype(k) ,'   Reservoir Type'
          do kk=1,npts(k)
            write(6,'(i5,10f8.1)') kk,
     &               rvelev(k,kk),rvarea(k,kk),rvstor(k,kk),
     &               rvout(k,kk,1), rvout(k,kk,2),
     &               rvout(k,kk,3), rvout(k,kk,4)
          end do
          write(6,'(f8.1,a)') slimit(k,4),'   Op Storage'
        end do
      end if

      ! read storm water diversion scheme.

      call OpnIF2(1,92, ifile,cstar)    
      
      ntransfer = 0
      pcap = 0.
      inquire(unit=92, opened = okfile)
      if(okfile) then
        write(6,'(/,a)') ' Water Transfer Data'
        call skipline(92,'*')
        Do it=1,10                                                     !maximum transfer
          Read(92,'(i4,4x,2a8,10f8.0)') itt, ch8(1), ch8(2),
     &                alpha(it), dperc(it), dmin(it), rv1
            ! transfer id number, 
            ! source and target recharge basin link name,
            ! fraction of target basin storage to fill with import
            ! target basin daily perc rate
            ! minimum source basin storage to transfer water
            ! transfer facility (pipe or canal) ID number
 
          if (itt.le.0) goto 72           ! need a blank line 
          if (itt.ne.it) then
            write(*,'(a)') ' Water transfer definition data error.'
            write(*,'(a)') ' Stop execution 1'
            stop
          end if
          ipipe(it) = rv1
          ntransfer = it 
          Call ctrim8r(ch8(1))
          Call ctrim8r(ch8(2))
          ! find source reservoir and link 
          Do k=1,nres
            j = lkatrv(k)
            if(lkname(j).eq.ch8(1)) then
              lkatts(it) = j
              kratts(it) = k ! rvname(kratet(it)) is reservoir name
              goto 70 
            end if
          end do
          write(*,'(a)') ' Water transfer definition data error.'
          write(*,'(i4,2a8)') itt,ch8(1),ch8(2)
          write(*,'(a)') ' Stop execution 2'
          stop
 70       continue
          Do k=1,nres
            j = lkatrv(k)
            if(lkname(j).eq.ch8(2)) then
              lkatte(it) = j
              kratte(it) = k
              goto 71
            end if
          end do
          write(*,'(a)') ' Water transfer definition data error.'
          write(*,'(i4,2a8)') itt,ch8(1),ch8(2)
          write(*,'(a)') ' Stop execution 3'
          stop
 71       continue
          ! write output
          write(6,'(i4,4x,2a8,3f8.2,I8)')  it,ch8(1),ch8(2),
     &                  alpha(it), dperc(it), dmin(it), ipipe(it)
          write(6,'(8X,2a32)')  rvname(kratts(it)), rvname(kratte(it))
        end do  ! it
     
        ! read transfer pipe data
        write(6,'(/,a)') ' Input data for water transfer facilities'
 72     do ip=1,10
          read(92,'(i4,4x,10f8.0)') kp, rtemp(1), rtemp(2)
          if(kp.le.0) goto 73
          pcap(kp) = rtemp(1)
          pmcap(kp) = rtemp(2)
          ! write pipe information to output file
          write(6,'(i8,2f8.2)') kp,pcap(kp), pmcap(kp)
          pcap(kp)  = pcap(kp)  * 1.98347
          pmcap(kp) = pmcap(kp) * 1.98347
        end do
        do it=1,ntransfer
          ip = ipipe(it)
          if(pcap(ip).le.0.) then
            write(*,'(a)') ' Water transfer definition data error.'
            write(*,'(2i8)') it,ip
            write(*,'(a)') ' Stop execution 4'
            stop
          end if
        end do
 73     close(92)

        read(1,'(a)') ffile
        if(ffile.eq.bl50) goto 74
        open(92,file=ffile,status='UNKNOWN')
        write(92,'(a,a)') 'Main input data file      = ',file0
        write(92,'(a)') 'Water Transfer Output File'
        write(92,'(a,a)') 'Water Transfer Definition = ',ffile
        write(92,'(/,a8,10i8,/)')      'YearMnDy',(it,it=1,ntransfer)
 74     continue        
      else 
        write(6,'(a)') 'No transfer of water between recharge basins'
      end if ! OKfile         
      
c
c     start routing by calling each link one at a time and storing
c     one time period at a time
c
      t=0
 1006 format(i4,2i2,25f8.0,/, (8x, 25f8.0))

c
c     ============================= Time Loop =============================
c
 100  t=t+1

      qout = 0.             ! (MXLK,0:MXPS)     flow out of nodes
      qlout = 0.            ! (MXLK,0:MXPS)     flow at end of links
      cqout = 0.            ! (MXLK,4,0:MXPS)   concentration out of nodes
      cqlout = 0.           ! (MXLK,4,0:MXPS)   concentration at end of links
      qloss = 0.            ! (MXLK,0:MXPS)     percolation from links
      qevp = 0.             ! (MXLK)            evaporation from links
	  pdivf = 0.
	  adivf = 0.
c
      read(11,1006,end=999) iyr,imo,idy,(qha(k),k=1,nhat)
      if(imo.eq.ismo .and. idy.eq.1) then 
        write(*,'(i5,2i2)') iyr
      end if
      ieyr = iyr
        
      mday = mndays(iyr,imo)


      ! decide reservoir operation mode for the day
      call rvopmd(iyr,imo,idy,7)

      if(qopt) then
        ! note that only option 1 is allowed, otherwise modify this routine
        read(21,1006,end=999) jyr,jmo,jdy,(cqha(k,1),k=1,nhat)
        call cmpdate(iyr,imo,idy, 21,jyr,jmo,jdy)

        read(31,1006,end=999) jyr,jmo,jdy,(cqha(k,2),k=1,nhat)
        call cmpdate(iyr,imo,idy, 31,jyr,jmo,jdy)
      end if
      if(ldebug) write(*,'(a)') ' pass q0a'

      if(nbf.gt.0) then
        read(12,'(i4,2i2,30f8.0)',end=999)
     &                jyr,jmo,jdy,(qha(k),k=nhat+1,nhat1)
        call cmpdate(iyr,imo,idy, 12,jyr,jmo,jdy)
        if(ldebug) write(*,'(a)') ' pass q0a1'
        if(qopt) then
          ! assign TDS                      ! boundary inflow
          if(mopt(2,1) .eq.1) then
            read(21,'(i4,2i2,30f8.0)') jr,jmo,jdy,
     &                     (cqha(k,1),k=nhat+1,nhat1)
            call cmpdate(iyr,imo,idy, 21,jyr,jmo,jdy)     
          elseif(mopt(2,1) .eq. 2) then
            do k=nhat+1,nhat1
              cqha(k,1) = cfixed(k,1)
            end do
          else
            do k=nhat+1,nhat1
              cqha(k,1) = 0.
              if(qha(k).gt.0.) then
                ! TDS
                flow = qha(k) * flowconv(k,1)
				ilu=ptrc(k,1)
                nt = ntconc(ilu,1)                                    
                do i=1,nt                                             
                  if(flow .le. tconc(ilu,1,i,1)) then            
                    conc = tconc(ilu,1,i-1,2)                       
     &                   + (tconc(ilu,1,i,2) - tconc(ilu,1,i-1,2))
     &                   / (tconc(ilu,1,i,1) - tconc(ilu,1,i-1,1))
     &                   * (flow             - tconc(ilu,1,i-1,1))
                    goto 31                                           
                  end if                                              
                end do                                                
                write(*,'(a)') ' Error in conc interpolation'         
                write(*,'(3i5,e12.4)') iha,ilu,1,flow
                stop  
 31             continue                                                
                cqha(k,1) = conc
                if(ldebug) write(*,'(a)') ' pass q0a2'
              end if
            end do
          end if ! mopt(2,1)
          ! assign TIN
          if(mopt(2,2) .eq.1) then                ! boundary inflow
            read(31,'(i4,2i2,10f8.0)') jr,jmo,jdy,
     &                     (cqha(k,2),k=nhat+1,nhat1)
            call cmpdate(iyr,imo,idy, 31,jyr,jmo,jdy)     
          elseif(mopt(2,2) .eq. 2) then
            do k=nhat+1,nhat1
              cqha(k,2) = cfixed(k,2)
            end do
          else
            do k=nhat+1,nhat1
              cqha(k,2) = 0.
              if(qha(k).gt.0.) then
                ! TIN
                flow = qha(k) * flowconv(k,2)
				ilu=ptrc(k,2)
                nt = ntconc(ilu,2)                                    
                do i=1,nt                                             
                  if(flow .le. tconc(ilu,2,i,1)) then            
                    conc = tconc(ilu,2,i-1,2)                       
     &                   + (tconc(ilu,2,i,2) - tconc(ilu,2,i-1,2))
     &                   / (tconc(ilu,2,i,1) - tconc(ilu,2,i-1,1))
     &                   * (flow             - tconc(ilu,2,i-1,1))
                    goto 32                                           
                  end if                                              
                end do                                                
                write(*,'(a)') ' Error in conc interpolation'         
                write(*,'(3i5,e12.4)') iha,ilu,2,flow
                stop  
 32             continue                                                
                cqha(k,2) = conc
                if(ldebug) write(*,'(a)') ' pass q0a2'
              end if
            end do
          end if ! mopt(2,2)
        end if ! qopt
      end if ! nbf
      if(ldebug) write(*,'(a)') ' pass q1a'

      ! remove lid flow
      do k=1,nhat
        flwlid(k) = min(qha(k),flid(k))
        qha(k) = qha(k) - flwlid(k)
      end do
      if(lidfile.ne.bl50) write(14,
     &        '(i4,2i2.2,25f8.2,/,(8x,25f8.2))')
     &          iyr,imo,idy, (flwlid(iha),iha=1,nhat)


      if(nps.gt.0) then
        read(13,'(i4,2i2,50f8.0)',end=999)
     &                jyr,jmo,jdy,(qha(k),k=nhat1+1,nhat2)
        call cmpdate(iyr,imo,idy, 13,jyr,jmo,jdy)
        if(qopt) then
          ! TDS
          if(mopt(3,1).eq.1) then
            read(23,'(i4,2i2,50f8.0)',end=999)
     &                jyr,jmo,jdy,(cqha(k,1),k=nhat1+1,nhat2)
            call cmpdate(iyr,imo,idy, 23,jyr,jmo,jdy)
          else
            do k=nhat1+1,nhat2
              cqha(k,1) = cfixed(k,1)
            end do
          end if
          ! TIN
          if(mopt(3,2).eq.1) then
            read(33,'(i4,2i2,50f8.0)',end=999)
     &                jyr,jmo,jdy,(cqha(k,2),k=nhat1+1,nhat2)
            call cmpdate(iyr,imo,idy, 33,jyr,jmo,jdy)
          else
            do k=nhat1+1,nhat2
              cqha(k,2) = cfixed(k,2)
            end do
          end if
        end if
      end if

	  if(npdd.gt.0) then
	    read(15,'(i4,2i2,30f8.0)') jyr,jmo,jdy,(pdivf(k),k=1,npdd)
		call cmpdate(iyr,imo,idy, 15,jyr,jmo,jdy)
	  end if
	  
      if(ldebug) write(*,'(a)') ' pass q1b'

      read(5,'(i4,2i2,10f8.0)',end=999) jyr,jmo,jdy,
     &                                  (evapor(k),k=1,newts)
      call cmpdate(iyr,imo,idy,  5,jyr,jmo,jdy)


      do 500 j=1,nlinks ! -----------------------------------------

      if(convtype(j).eq.5) then
c
c       dummy split linkc       do nothing
c
        goto 400
      end if

      iun = usn(j)

      ! historical or planned diversion as negative inflow
      rval3 = 0.

      ! add point and non-point source flows
      if (nha(iun).gt.0) then
        do k=1,nha(iun)
          kha = nodha(iun,k)
          kps = kha - nhat1
          rval0 = qha(kha)*rmult(iun,k)
          rval1 = rval0 * cqha(kha,1)
          rval2 = rval0 * cqha(kha,2)
          if(rval0.gt.0.) then
            qout(iun,0) = qout(iun,0) + rval0
            if(qopt) cqout(iun,3,0) = cqout(iun,3,0) + rval1        ! total mass
            if(qopt) cqout(iun,4,0) = cqout(iun,4,0) + rval2
            if(kps.gt.0) then
              qout(iun,kps) = qout(iun,kps) + rval0
              if(qopt) cqout(iun,3,kps) = cqout(iun,3,kps) + rval1  ! mass by kps
              if(qopt) cqout(iun,4,kps) = cqout(iun,4,kps) + rval2
            end if
          else
            ! historical or planned diversion
            rval3 = rval3 - rval0
          end if
        end do
      end if
      if(ldebug) write(*,'(a)') ' pass q2a'

      ! add upstream link inflows
      if(nincl(iun).gt.0) then
        do kk=1,nincl(iun)
          ilk = incoml(iun,kk)
          do kps=0,nps
            qout(iun,kps) = qout(iun,kps) + qlout(ilk,kps)
            if(qopt) cqout(iun,3,kps) = cqout(iun,3,kps)
     &                           + qlout(ilk,kps)*cqlout(ilk,1,kps)
            if(qopt) cqout(iun,4,kps) = cqout(iun,4,kps)
     &                           + qlout(ilk,kps)*cqlout(ilk,2,kps)
          end do
        end do
      end if
      if(ldebug) write(*,'(a)') ' pass q2b'

      ! calculate concentration
      if(qopt) then
        do kps=0,nps
          if(qout(iun,kps).gt.0.) then
            cqout(iun,1,kps) = cqout(iun,3,kps) / qout(iun,kps)
            cqout(iun,2,kps) = cqout(iun,4,kps) / qout(iun,kps)
          end if
        end do
      end if
      if(ldebug) write(*,'(a)') ' pass q2c'

      ! subtract any historical or planned diversion
c     if(rval3.gt.0.) then
      if(rval3.gt.0.1 ) then            !4/17/07
        rval0 = qout(iun,0)
        qout(iun,0) = max(rval0-rval3,0.)
c            cqout(iun,1,0) = 
c            cqout(iun,2,0) = 
        rval1 = qout(iun,0) / rval0
        do kps=1,nps
          qout(iun,kps) = qout(iun,kps) * rval1
c            cqout(iun,1,kps)
c            cqout(iun,2,kps)
         end do
      end if

      if(ldebug) write(*,'(a,2i5)') '  Pass q2d',j,convtype(j)

      if(ntransfer.gt.0) then   
        if(ldebug) write(*,'(a,2i5)') '  Pass q2e'
        do it = 1,ntransfer
          jun = usn(lkatte(it))
          if(jun.eq.iun) then
            qout(iun,0) = qout(iun,0) + pimp(it)/1.98347
c            cqout(iun,1,0) = 
c            cqout(iun,1,0) = 
          end if
        end do
      end if
      
      goto (150,155,200,130,400,140) convtype(j)  ! ===================

 130  continue
      if(ldebug) write(*,'(a)') '  Pass q3a'
c     diversion link
      idiv = ptrdv(j)
      do kk=2,ndpts(idiv)
        if(qout(iun,0).le.dflow(idiv,kk)) go to 490
      end do
      kk=ndpts(idiv)
 490  slopei = (qout(iun,0)-dflow(idiv,kk-1))
     &        /(dflow(idiv,kk)-dflow(idiv,kk-1))
      div1= divout(idiv,kk-1)
     &    + slopei*(divout(idiv,kk)-divout(idiv,kk-1))
      div1 = max(0.,div1)
      div2 = max(0., qout(iun,0)-div1)
c            cqout(iun,1,0) = 
c            cqout(iun,1,0) = 
 

      r1 = 0.
      if(qout(iun,0).gt.0.) r1 = div1 / qout(iun,0)
      r2 = 1. - r1

c     send the flow to downstream links
      ilk1 = linkdv(idiv,1)          ! diversion
      ilk2 = linkdv(idiv,2)          ! main stream
      do kps=0,nps
        if(ilk1.gt.0) then
          qlout(ilk1,kps) = qout(iun,kps) * r1
          if(qopt) cqlout(ilk1,1,kps) = cqout(iun,1,kps)
          if(qopt) cqlout(ilk1,2,kps) = cqout(iun,2,kps)
        end if
        if(ilk2.gt.0) then
          qlout(ilk2,kps) = qout(iun,kps) * r2
          if(qopt) cqlout(ilk2,1,kps) = cqout(iun,1,kps)
          if(qopt) cqlout(ilk2,2,kps) = cqout(iun,2,kps)
        end if
      end do
      if(divfile.ne.bl50 .and. qout(iun,0).gt.0.)
     &     write(10,'(6i5,3f10.1)') imo,idy,jyr,
     &                        j,ilk1,ilk2,
     &                        qout(iun,0),div1,div2
      if(ldebug) write(*,'(a)') '  Pass q3b'
      goto 400


      ! convtype 6 is a special case of convtype 1
      ! flow-width rating table is used for this type
 140  continue
      ich = ptrch(j)
      flow = qout(iun,0)
      if(ldebug) write(*,'(a,2i5,f10.2))') '  Pass2',j,ich,flow
c      width(j) = rintp(0,ich,flow)
      width(j) = rintpm(0,ich,flow)
      qloss(j,0) = min(qout(iun,0), (length(j)*width(j)*dpb(j)) )
      if(ldebug) write(*,'(a,5f10.2))') '  Pass2',width(j),qloss(j,0)
      go to 165

c     check to see if it is necessary to compute stage
c     test is different for storage and conveyance links
c        based on non-postive slope for conveyance only link
c        always compute stage if storage link

 150  continue
      if(qout(iun,0).le.0.00001) goto 160
      if(slop(j).le.0) go to 160 ! this is lined channel,  just move water
      call stagec(slop(j),mannn(j),b(j),z(j),dpb(j),dpz(j),
     1  qout(iun,0),depth,unitperc)

      if(depth.lt.0.) then
         write(6,'(i5,10f12.4)') j,slop(j),mannn(j),b(j),z(j),dpb(j),
     &                           unitperc
        stop
      end if
      qloss(j,0) = min(qout(iun,0),unitperc*length(j))
      go to 165

 155  call stagep()
 160  qloss(j,0) = 0.
 165  r1 = 0.
      if(qout(iun,0).gt.0.) r1 = qloss(j,0) / qout(iun,0)

      do kps=0,nps
        qloss(j,kps) = qout(iun,kps) * r1
        qlout(j,kps) = qout(iun,kps) - qloss(j,kps)
        if(qopt) cqlout(j,1,kps) = cqout(iun,1,kps)
        if(qopt) cqlout(j,2,kps) = cqout(iun,2,kps) * exp(-1.*rxn(j))
      end do
      if(ldebug) write(*,'(a)') '  Pass q4a'
      go to 400

c
c     reservoir routing using rating curves
c
 200  continue
      if(ldebug) write(*,'(a)') '  Pass q5a'
      do k=1,nres
        if(lkatrv(k).eq.j) go to 210
      end do
      go to 990
 210  continue

      ! reset variables
      ilk1 = linkrv(k,1)
      ilk2 = linkrv(k,2)
      evapl(k) = 0.
      percl(k) = 0.
      out1l(k) = 0.
      out2l(k) = 0.

      qin(k) = qout(iun,0)
      evapl(k) = 0.
      do kk=1,newts
        evapl(k) = evapl(k) + evwts(k,kk)*evapor(kk)/12.
      end do
c

      qin_af = cf1*qin(k)
                                    ! storl(k,0) is previous day total storage
      stor0 = storl(k,0) + qin_af                        ! save for concentration calculation
      cstor1 = 0.
      cstor2 = 0.
      if(qopt .and. stor0.gt.0.) then
        cstor1 = (storl(k,0) * cstorl(k,1,0)             ! initial concentration
     &            + qin_af * cqout(iun,1,0)) / stor0
        cstor2 = (storl(k,0) * cstorl(k,2,0)
     &            + qin_af * cqout(iun,2,0)) / stor0
      end if



c     call rv_sim(k,qin_af)                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call rvsim2(k,qin_af)                ! rvsim2 has storm mode operation
                                    ! after rvsim, storl(k,0) is today's end storage
      r1 = 0.
      r2 = 0.
      r3 = 0.
      r4 = 0.
      r5 = 0.
      r6 = 0.
      r7 = 0.
	  
	  ! planned diversion                                   !v6a2
      if(npdd.gt.0) then 
	    kpd=kpdiv(j)                                        !v6a2
	    if(kpd.gt.0.) then                                  !v6a2
	      pdivaf = pdivf(kpd)*cf1                           !v6a2
	      if(pdivaf.gt.storl(k,0)) pdivaf = storl(k,0)      !v6a2
	   	  adivf(kpd) = pdivaf/cf1                           !v6a2
	   	  storl(k,0) = storl(k,0) - pdivaf                  !v6a2
        end if                                              !v6a2 
	  end if
      if(qopt .and. stor0.gt.evapl(k)) then
        r1 = (stor0 - evapl(k)) / stor0
        r2 = 1. / r1
        if(rxn(j).ge.0.) then
          r3 = r2 * exp(-rxn(j))
        else
          if(cstor2.gt.0.) r3 = -rxn(j) / cstor2
        end if
        !!! r4 = percl(k) / stor0
        r5 = out1l(k) / stor0
        r6 = out2l(k) / stor0
        r7 = storl(k,0) / stor0
      end if
      if(qopt) then
        if(stor0 .gt.0.) then
          cstorl(k,1,0) = cstor1 * r2
          cstorl(k,2,0) = cstor2 * r3
        else
          cstorl(k,1,0) = 0.
          cstorl(k,2,0) = 0.
        end if

      end if ! qopt

      if(ldebug) write(*,'(a)') ' pass q5b'

c     send the flow to downstream links
      if(ilk1.gt.0) then
        qlout(ilk1,0) = out1l(k)
        if(qopt) cqlout(ilk1,1,0) = cstorl(k,1,0)
        if(qopt) cqlout(ilk1,2,0) = cstorl(k,2,0)
      end if
      if(ilk2.gt.0) then
        qlout(ilk2,0) = out2l(k)
        if(qopt) cqlout(ilk2,1,0) = cstorl(k,1,0)
        if(qopt) cqlout(ilk2,2,0) = cstorl(k,2,0)
      end if

c     write(*,'(3i6,f10.1,f10.2)')
c    &      j,k,ilk1,cqlout(ilk1,1,0),cqlout(ilk1,2,0)  ! delete cqlout
c     write(*,'(3i6,f10.1,f10.2)')
c    &      j,k,ilk2,cqlout(ilk2,1,0),cqlout(ilk2,2,0)  ! delete cqlout

      ! routine for point source components
      do kps=1,nps
        qin_af = cf1 * qout(iun,kps)
        stor0 = storl(k,kps) + qin_af
        cstor1 = 0.
        cstor2 = 0.
        if(qopt .and. stor0.gt.0.) then
          cstor1 = (storl(k,kps) * cstorl(k,1,kps)
     &              + qin_af * cqout(iun,1,kps)) / stor0
          cstor2 = (storl(k,kps) * cstorl(k,2,kps)
     &              + qin_af * cqout(iun,2,kps)) / stor0
        end if
        !!! percl(k,kps) = stor0 * r4
        storl(k,kps) = stor0 * r7
        if(qopt) cstorl(k,1,kps) = cstor1 * r2
        if(qopt) cstorl(k,2,kps) = cstor2 * r3
        if(ilk1.gt.0) then
          qlout(ilk1,kps) = stor0 * r5
          if(qopt) cqlout(ilk1,1,kps) = cstorl(k,1,kps)
          if(qopt) cqlout(ilk1,2,kps) = cstorl(k,2,kps)
        end if
        if(ilk2.gt.0) then
          qlout(ilk2,kps) = stor0 * r6
          if(qopt) cqlout(ilk2,1,kps) = cstorl(k,1,kps)
          if(qopt) cqlout(ilk2,2,kps) = cstorl(k,2,kps)
        end if
        if(qopt .and. stor0.gt.0.) then
          cstorl(k,1,kps) = (storl(k,kps) * cstorl(k,1,kps)
     &              + qin_af * cqout(iun,1,kps) ) / stor0 * r2
          cstorl(k,2,kps) = (storl(k,kps) * cstorl(k,2,kps)
     &              + qin_af * cqout(iun,2,kps) ) / stor0 * r3
        else
          cstorl(k,1,kps) = 0.
          cstorl(k,2,kps) = 0.
        end if
        if(ldebug) write(*,'(a)') ' pass q5c'

        if(ilk1.gt.0) then
          qlout(ilk1,kps) = stor0 * r5
          if(qopt) cqlout(ilk1,1,kps) = cstorl(k,1,kps)
          if(qopt) cqlout(ilk1,2,kps) = cstorl(k,2,kps)
        end if
        if(ilk2.gt.0) then
          qlout(ilk2,kps) = stor0 * r6
          if(qopt) cqlout(ilk2,1,kps) = cstorl(k,1,kps)
          if(qopt) cqlout(ilk2,2,kps) = cstorl(k,2,kps)
        end if
        !!! percl(k,kps) = stor0 * r4
        storl(k,kps) = stor0 * r7
      end do

      qloss(j,0) = percl(k) / cf1
      qevp(j) = evapl(k)

 400  continue    ! continuation after each link
 500  continue    ! j=1,nlinks loop

 2000 format(5i5,7f10.1)
 2090 format(3i5,30f8.1 / 3(15x,30f8.1 / ))
 2091 format(3i5,30i8 / 3(15x,30i8 / ))


      do j = 1,nlinks
        if(qlmax(j).lt.qlout(j,0)) then
          jqlmax(j) = jyr*10000 + imo*100 + idy
          qlmax(j) = qlout(j,0)
        end if
      end do

      if(ldebug) write(*,'(a)') '  Pass q6a' !delete ---------------------------
      if(nplk.gt.0) then
        do i=1,nplk
          if(iplk(i,2).eq.1) then
            rtemp(i) = qlout(iplk(i,1),0)
          else
            rtemp(i) = qout(iplk(i,1),0)
          end if
        end do
        write(61,'(i4.4,2i2.2,30f8.1)') iyr,imo,idy,
     &                            (rtemp(i),i=1,nplk)
     
        if(qopt) then
          do i=1,nplk
            if(iplk(i,2).eq.1) then
              rtemp(i) = cqlout(iplk(i,1),1,0)
              rtemp2(i) = cqlout(iplk(i,1),2,0)
            else
              rtemp(i) = cqout(iplk(i,1),1,0)
              rtemp2(i) = cqout(iplk(i,1),2,0)
            end if
          end do
          write(62,'(i4.4,2i2.2,30f8.2)') iyr,imo,idy,
     &                            (rtemp(i),i=1,nplk)
          write(63,'(i4.4,2i2.2,30f8.2)') iyr,imo,idy,
     &                            (rtemp2(i),i=1,nplk)
        end if
      end if !nplk
      if(ldebug) write(*,'(a)') '  Pass q6b' !delete ---------------------------

      if(nqls.gt.0) then
        do i=1,nqls
c          xqloss(i) = xqloss(i) + qloss(ipql(i),0) * cf1
          if(ipql(i) .gt. 0) then             ! 1/29/19
c            xqloss(i) = xqloss(i) + qloss(kpql(i),0)                    ! 1/29/19
c            xqevp(i) = xqevp(i) + qevp(kpql(i))                         ! 1/29/19
            xqloss(i) = xqloss(i) + qloss(ipql(i),0)                     ! 1/29/19
            xqevp(i) = xqevp(i) + qevp(ipql(i))			                 ! 1/29/19
            if(qopt) then
c              cxqloss(i,1) = cxqloss(i,1) +                              ! 1/29/19
c     &                       qloss(kpql(i),0) * cstorl(kpql(i),1,0)      ! 1/29/19
c              cxqloss(i,2) = cxqloss(i,2) +                              ! 1/29/19
c     &                       qloss(kpql(i),0) * cstorl(kpql(i),2,0)      ! 1/29/19
               cxqloss(i,1) = cxqloss(i,1) +                              ! 1/29/19
     &                       qloss(ipql(i),0) * cstorl(kpql(i),1,0)       ! 1/29/19
              cxqloss(i,2) = cxqloss(i,2) +                               ! 1/29/19
     &                       qloss(ipql(i),0) * cstorl(kpql(i),2,0)       ! 1/29/19
            end if 
		  else               ! 1/29/19      
			xqloss(i) = 0.   ! 1/29/19
		  end if             ! 1/29/19
        end do
        if(idy.eq.mday) then
          do i=1,nqls
            if(qopt) then
              if(xqloss(i).gt.0.) then
                cxqloss(i,1) = cxqloss(i,1) / xqloss(i)
                cxqloss(i,2) = cxqloss(i,2) / xqloss(i)
              else
                cxqloss(i,1) = 0.
                cxqloss(i,2) = 0.
              end if
            end if    
            xqloss(i) = xqloss(i) * cf1
          end do
          write(81,'(i4,i2.2,(50f8.2))') iyr,imo,
     &                              (xqloss(i),i=1,nqls)

          if(qopt) then
            write(82,'(i4,i2.2,(50f8.1))') iyr,imo,            
     &                              (cxqloss(i,1),i=1,nqls)
            write(83,'(i4,i2.2,(50f8.2))') iyr,imo,            
     &                              (cxqloss(i,2),i=1,nqls)
          end if
          write(84,'(i4,i2.2,(50f8.2))') iyr,imo, 
     &                              (xqevp(i),i=1,nqls)     

          do i=1,nqls
            xqloss(i) = 0.
            cxqloss(i,1) = 0.
            cxqloss(i,2) = 0.
            xqevp(i) = 0.
          end do
        end if ! idy
        
        inquire(unit=85, opened = okfile)
        if(okfile) write(85,'(i4,2i2.2,(50f8.2))') iyr,imo,idy,
     &                  (storl(kpql(i),0),i=1,nqls)
      end if !nqls

      if(npdd.gt.0) then 
        inquire(unit=91, opened = okfile)
        if(okfile) write(91,'(i4,2i2.2,(50f8.2))') iyr,imo,idy,
     &                    (adivf(kpd),kpd=1,npdd)
      end if
      if(ldebug) write(*,'(a)') '  Pass q6c' !delete ---------------------------

      if(nprch.gt.0) then
        do j=1,nlinks
          ii = linkrch(j)
          if(ii.gt.0 .and. ii.le.30) then
            iun = usn(j)
c            do kps=0,nps
c              archp(ii,0,kps,iyr,imo) = archp(ii,0,kps,iyr,imo)
c     &                        + qloss(j,kps)
c              if(qopt) archp(ii,1,kps,iyr,imo) = archp(ii,1,kps,iyr,imo)
c     &                        + qloss(j,kps) * cqlout(j,1,kps)
c              if(qopt) archp(ii,2,kps,iyr,imo) = archp(ii,2,kps,iyr,imo)
c     &                        + qloss(j,kps)
c     &                        * (cqlout(j,2,kps) + cqout(iun,2,kps))/2.
c            end do
              archp(ii,0,iyr,imo) = archp(ii,0,iyr,imo)
     &                        + qloss(j,0)
              if(qopt) archp(ii,1,iyr,imo) = archp(ii,1,iyr,imo)
     &                        + qloss(j,0) * cqlout(j,1,0)
              if(qopt) archp(ii,2,iyr,imo) = archp(ii,2,iyr,imo)
     &                        + qloss(j,0)
     &                        * (cqlout(j,2,0) + cqout(iun,2,0))/2.
          end if
        end do
      end if
      if(ldebug) write(*,'(a)') '  Pass q6d' !delete ---------------------------

      if(npnd.gt.0) then
        do i=1,npnd
          ii = ipnd(i)
          aflow(i,0,iyr,imo) = aflow(i,0,iyr,imo)
     &                         + qout(ii,0)
          if(qopt) aflow(i,1,iyr,imo) = aflow(i,1,iyr,imo)
     &                         + qout(ii,0) * cqout(ii,1,0)
          if(qopt) aflow(i,2,iyr,imo) = aflow(i,2,iyr,imo)
     &                         + qout(ii,0) * cqout(ii,2,0)
        end do
      end if !npnd
      if(ldebug) write(*,'(a)') '  Pass 6e' !delete ---------------------------

      
      !at end of day, after routing all reaches is done, calculate possible transfer for following day
      if(ntransfer.gt.0) then         
        if(ldebug) write(*,'(a)') '  Pass 7a'   
        ttrans = 0.                   
        ptrans = 0.                   
        prcap = pcap   ! reset for all pipe capacity
        
        if(ldebug) write(*,'(a)') '  Pass 7b'                                      
        do it=1,ntransfer             
          kts = kratts(it)  ! reservoir serial number for transfer source
          kte = kratte(it)  !                         for transfer target
c          jts = lkatts(it)  ! link serial number for transfer source
c          jte = lkatte(it)  !                    for transfer target
c          iune = usn(jte)   ! upstream node for target reservoir
                                      
          ip = ipipe(it)    ! tranfer facility, pipe or canal or any transfer mechanism, ID number
          sttgt = storl(kte,0)     ! day-end storage of target reservoir
          stmax = slimit(kte,4)    ! maximum storage of target reservoir
          stsrc = storl(kts,0)        
          ! dmin(it) minimum storage at source to transfer water
          rv1 = max(stmax*alpha(it)-sttgt,0.)+dperc(it)    ! target reservoir condition
          rv1 = min(rv1,prcap(ip))                  ! pipe condition
          rv1 = min(rv1,max(stsrc-dmin(it),0.))            ! source reservoir condition
          ! now rv1 is daily transfer for next day
          ! now adjust storage and capacity
          ttrans(it) = rv1            
          prcap(ip) = prcap(ip) - rv1                      ! reduce remaining pipe capacity
          pimp(it) = rv1                                   ! save to load at USN        
          storl(kts,0) = storl(kts,0) - rv1                ! reduce source storage 
          ptrans(ip) = ptrans(ip) + ttrans(it)
        end do
        if(ldebug) write(*,'(a)') '  Pass 7e'        
        inquire(unit=92, opened = okfile)
        if(okfile) write(92,'(i4,2i2.2,10f8.2)' ) iyr,imo,idy,
     &                                 (pimp(it),it=1,ntransfer)
      end if                          
      
      
      
      
      go to 100

 990  write(*,9900)
 9900 format('program stopped because link ', i5,' reservoir data'
     1 ,' is missing')
      stop

 999  write(*,'(a,i4,2i2.2)') ' Stop on ',iyr,imo, idy

      write(6,'(/,a,/)') ' Maximum Flow at Each Link'
      do j=1,nlinks
        write(6,'(i5,2x,a8,i10,f10.2)') j,lkname(j),jqlmax(j),qlmax(j)
      end do
c      close(2)
c      close(3)
c      close(4)
c      close(6)
      close(11)
      close(12)
      close(13)
c      close(21)
c      close(22)
c      close(23)
c      close(31)
c      close(32)
c      close(33)
c      close(61)
c      close(62)
c      close(63)

      if(npnd.gt.0) then !----------------------------------------------
c     do iyr=1920,2011
	  do iyr=isyr,ieyr
        do imo=1,12
          do ii=1,npnd
            if(aflow(ii,0,iyr,imo).le.0.) aflow(ii,0,iyr,imo)=0.
            if(aflow(ii,0,iyr,imo).gt.0.049) then             !!!!!!!!!!!!!!!!!0.) then
               aflow(ii,1,iyr,imo) = aflow(ii,1,iyr,imo)/
     &                               aflow(ii,0,iyr,imo)     ! TDS mg/L
               aflow(ii,2,iyr,imo) = aflow(ii,2,iyr,imo)/
     &                               aflow(ii,0,iyr,imo)     ! TIN mg/L
            else
               aflow(ii,1,iyr,imo) = 0.
               aflow(ii,2,iyr,imo) = 0.
            end if                 
            aflow(ii,0,iyr,imo) = aflow(ii,0,iyr,imo) * cf1     !(af)
c            aflow(ii,1,iyr,imo) = aflow(ii,1,iyr,imo) * cf2    !tons
c            aflow(ii,2,iyr,imo) = aflow(ii,2,iyr,imo) * cf2    !tons
          end do
        end do
      end do

c       write(64,'(a)') ' Total monthly flow passing the nodes'
c       write(64,'(a)') ' Units (ac-ft)'

c       do ify=1921,2011
        do ify=isyr,ieyr
          do ifm=1,12
            if(ifm.le.6) then
              icy = ify - 1
              icm = ifm+6
            else
              icy = ify
              icm = ifm-6
            end if
            write(64,'(i4,i2,30f8.1)')
     &           icy,icm,(aflow(ii,0,icy,icm),ii=1,npnd)
            if(qopt) then
              write(65,'(i4,i2,30f8.1)')
     &             icy,icm,(aflow(ii,1,icy,icm),ii=1,npnd)
              write(66,'(i4,i2,30f8.2)')
     &             icy,icm,(aflow(ii,2,icy,icm),ii=1,npnd)
            end if
          end do
        end do
      end if

c     end if !(npnd.gt.0) ----------------------------------------------

      if(nqls.gt.0) close(81)
      if(nqls.gt.0) close(82)
      if(nqls.gt.0) close(85)

      if(nprch.gt.0) then !---------------------------------------------
c       do iyr=1920,2011                                                 !7/8/05
        do iyr=isyr,ieyr                                                 !ro_v6a1
          do imo=1,12                                                    !7/8/05
c            do kps=0,nps                                                 !7/8/05
c              do ii=1,nprch                                              !7/8/05
c                if(archp(ii,0,kps,iyr,imo).gt.0.) then
c                   archp(ii,1,kps,iyr,imo) = archp(ii,1,kps,iyr,imo) /        !(mg/L)
c     &                                       archp(ii,0,kps,iyr,imo)
c                   archp(ii,2,kps,iyr,imo) = archp(ii,2,kps,iyr,imo) /        !(mg/L)
c     &                                       archp(ii,0,kps,iyr,imo)
c                 else
c                   archp(ii,1,kps,iyr,imo) = 0.
c                   archp(ii,2,kps,iyr,imo) = 0.
c                 end if
c                 archp(ii,0,kps,iyr,imo) = archp(ii,0,kps,iyr,imo) * cf1     !7/8/05
cc                archp(ii,1,kps,iyr,imo) = archp(ii,1,kps,iyr,imo) * cf2  !7/8/05 (tons)
cc                archp(ii,2,kps,iyr,imo) = archp(ii,2,kps,iyr,imo) * cf2  !7/8/05 (tons)
c              end do                                                     !7/8/05
c            end do                                                       !7/8/05

              do ii=1,nprch                                              !7/8/05
                if(archp(ii,0,iyr,imo).gt.0.049) then
                   archp(ii,1,iyr,imo) = archp(ii,1,iyr,imo) /        !(mg/L)
     &                                   archp(ii,0,iyr,imo)
                   archp(ii,2,iyr,imo) = archp(ii,2,iyr,imo) /        !(mg/L)
     &                                   archp(ii,0,iyr,imo)
                 else
                   archp(ii,1,iyr,imo) = 0.
                   archp(ii,2,iyr,imo) = 0.
                 end if
                 archp(ii,0,iyr,imo) = archp(ii,0,iyr,imo) * cf1     !7/8/05
c                archp(ii,1,iyr,imo) = archp(ii,1,iyr,imo) * cf2     !7/8/05 (tons)
c                archp(ii,2,iyr,imo) = archp(ii,2,iyr,imo) * cf2     !7/8/05 (tons)
              end do                                                     !7/8/05

          end do                                                         !7/8/05
        end do                                                           !7/8/05


cccccc        if(qopt) write(68,'(a)') ' Total monthly TDS mass percolation'
cccccc        if(qopt) write(68,'(a)') ' Units|(tons)'
cccccc        if(qopt) write(69,'(a)') ' Total monthly TIN mass percolation'
cccccc        if(qopt) write(69,'(a)') ' Units|(tons)'

        kps = 0
        do ify=isyr,ieyr                                                  !ro_v6a1
          do ifm=1,12
            if(ifm.le.6) then
              icy = ify - 1
              icm = ifm+6
            else
              icy = ify
              icm = ifm-6
            end if
            write(67,'(i4,i2,10f8.1)')
     &           icy,icm,(archp(ii,0,icy,icm),ii=1,nprch)
            if(qopt) then
              write(68,'(i4,i2,10f8.1)')
     &             icy,icm,(archp(ii,1,icy,icm),ii=1,nprch)
              write(69,'(i4,i2,10f8.2)')
     &             icy,icm,(archp(ii,2,icy,icm),ii=1,nprch)
            end if
          end do
        end do
      end if

      close(67)
      if(qopt) close(68)
      if(qopt) close(69)

      stop
      end




