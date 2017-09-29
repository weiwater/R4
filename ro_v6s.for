CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                          C
C  Rainfall-Runoff-Router-Rootzone Model                                   C
C  Copy Right (C) 1996 - 2017 Wildermuth Environmental, Inc.               C
C                                                                          C
C  This program is free software: you can redistribute it and/or modify    C
C  it under the terms of the GNU General Public License as published by    C
C  the Free Software Foundation, either version 3 of the License, or       C
C  (at your option) any later version.                                     C
C                                                                          C
C  This program is distributed in the hope that it will be useful,         C
C  but WITHOUT ANY WARRANTY; without even the implied warranty of          C
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           C
C  GNU General Public License for more details.                            C
C                                                                          C
C  You should have received a copy of the GNU General Public License       C
C  along with this program.  If not, see <http://www.gnu.org/licenses/>.   C
C                                                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C 10/20/2015 ro_v6s.for - cleaned version
C
C 7/21/2014 modified from ro_v6a.for to print daily storage in recharge basins
C

C Files opened
c  Unit  File description
c    1    main input file
c    2    link file
c    3    node file
c    4    open channel flow-width data (closed after data are read)
c         TDS table file (closed after data are read)
c         TIN table file (closed after data are read)
c    5    daily davporation data file
c    6    main output file
c    7    storm prediction rainfall data file                               !r5
c    9
c   10    diversion output file (optional)
c   11    runoff model output flow data file
c   12    boundary inflow data file
c   13    point source flow data file
c   21    runoff model output TDS data file
c   22    boundary inflow TDS data file
c   23    point source TDS data file
c   31    runoff model output TIN data file
c   32    boundary inflow TIN data file
c   33    point source TIN data file
c   61    output file to print flow at user specified link or node
c   62    output file to print TDS  at user specified link or node
c   63    output file to print TIN  at user specified link or node
c   64    output file to print flow at user spefied node (for storm water and each point source)
c   65    output file to print TDS  at user spefied node (for storm water and each point source)
c   66    output file to print TIN  at user spefied node (for storm water and each point source)
c   67    output file to print percolation in user spedifed reaches (for storm water and each point source)
c   68    output file to print TDS in percolation in user spedifed reaches (for storm water and each point source)
c   69    output file to print TIN in percolation in user spedifed reaches (for storm water and each point source)
c   71    output file to print percolation at user-specfied reach (older version, not used in current version)
c   72    WIDTH.OUT --- for type 6 conveyance link, print width               (commented out in current version)
c   73    QLOSS.OUT --- for type 6 conveyance link, print percolation (qloss) (commented out in current version)

      include 'RO_DIM_4.MAX'

      logical      okfile
      character*50 file0,lfile,nfile,hafile,bifile,psfile,
     &             efile,ofile,divfile,bl50,chfile,
     &             opfile(3,2),ffile,lidfile,czfile,
     &             prfile        ! precipitation file as storm forecaster
      character*30 fmtin

      real         rtemp(MXHA),rtemp2(MXHA)     ! temporary variables

c---- general model control properties
      include 'GenParm.VAR'

c---- link/node properties
      include 'LKprpt.var'

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

      ierror = .FALSE.
      jerror = 0
      bl50 = '                                                  '
      velocity = 1.0

c     initialize variables
      call vrinit()

      archp = 0.                     ! (MXRH,0:2,MXPS,1949:2006,12)
      aflow = 0.                     ! (MXRH,0:2,MXPS,1949:2006,12)
      kpql = 0    

c
c     read in control parameters from control file
c

      call getcon(file0,n1)

      open(1,file=file0,status='OLD')
      call skipline(1,'*')
      read(1,1003) newts,nhat,nbf,nps,idebug,iqopt
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
         ! ndffday, number of drydays to qualify as first flush
         ! rainss, minimum rainfall for significant storm


 1003 format(i8)


      call openifile(1, 2, lfile)         ! link definition file
      call openifile(1, 3, nfile)         ! node definition file
      call openifile(1, 4, ffile)         ! reach parameter file

      read(4,1003) ncnz
      do i=1,ncnz
        read(4,1003) j
      end do
      call skipline(4,'*')
      read(4,1003) nrainz
      do i=1,nrainz
        read(4,1003) j
      end do
      call skipline(4,'*')
      read(4,1003) nrchz
      do i=1,nrchz
        read(4,'(i8,4f8.0)') j,(rvrch(j,k),k=1,3)
        if(i.ne.j) then
          write(*,'(a)') ' Error in parameter file'
          stop
        end if
      end do
      close(4)


! flow input files
      call openifile(1,11,hafile)         ! runoff model output file
      call openifile(1,12,bifile)         ! recorded boundary inflow file
      call openifile(1,13,psfile)         ! recorded point source fileile
      call openifile(1,14,czfile)         ! ha area parameters
      read(14,'(a40)') fmtin
      do j=1,nhat
        read(14,fmtin,err=992,end=992) ii,ch8,
     &             areaha, i1,i2,r1,r2
        ! find HAID
      end do
 992  close(14)


      if(qopt) then
        ! TDS input files
        call openifile(1,21,opfile(1,1))
        if(opfile(1,1).ne.bl50) read(21,1003) mopt(1,1)
        ! only option type 1 is available

        call openifile(1,22,opfile(2,1))
        if(opfile(2,1).ne.bl50) read(22,1003) mopt(2,1)
        ! only option type 3 is available
        if(qopt .and. nbf.gt.0) then
          read(22,'(8x,30i8)') (ptrc(i,1),i=nhat+1,nhat1)
          read(22,'(8x,30f8.0)') (flowconv(i,1),i=nhat+1,nhat1)
        end if

        call openifile(1,23,opfile(3,1))
        if(opfile(2,1).ne.bl50) read(23,1003) mopt(3,1)
        ! only option types 1 and 2 are allowed
        if(mopt(3,1).eq.3) stop
        ! if option 1, do nothing
        ! if option 2, read constant concentration
        if(mopt(3,1).eq.2) read(23,'(8x,50f8.0)')
     &                              (cfixed(i,1),i=nhat1+1,nhat2)


        ! TIN input files
        call openifile(1,31,opfile(1,2))    ! TIN input data file
        if(opfile(1,2).ne.bl50) read(31,1003) mopt(1,2)
        ! only option type 1 is available

        call openifile(1,32,opfile(2,2))    ! TIN input data file
        if(opfile(2,1).ne.bl50) read(32,1003) mopt(2,2)
        ! only option type 3 is available
        if(qopt .and. nbf.gt.0) then
          if(mopt(2,2).ne.3) stop
          read(32,'(8x,30i8)') (ptrc(i,2),i=nhat+1,nhat1)
          read(32,'(8x,30f8.0)') (flowconv(i,2),i=nhat+1,nhat1)
        end if

        call openifile(1,33,opfile(3,2))    ! TIN input data file
        if(opfile(2,1).ne.bl50) read(33,1003) mopt(3,2)
        ! only option types 1 and 2 are allowed
        if(mopt(3,2).eq.3) stop
        ! if option 1, do nothing
        ! if option 2, read constant concentration
        if(mopt(3,2).eq.2) read(33,'(8x,50f8.0)')
     &                              (cfixed(i,2),i=nhat1+1,nhat2)
      end if


! channel rating curves

      call openifile(1, 4,chfile)         ! channel flow-width rating file
      nch = 0
      if(chfile.ne.bl50) then
        read(4,'(i8)') nch
        if(nch.gt.MXCH) then
          write(*,'(a)') ' Increase MXCH'
          stop
        end if
        do ich=1,nch
          read(4,'(a8,i6,50f6.0)') chname(ich),ntab(0,ich),
     &                             (table(0,ich,k,1),k=1,30)   ! flow
          read(4,'(8x,6x,50f6.0)') (table(0,ich,k,2),k=1,30)   ! width
          call ctrim8r(chname(ich))
        end do
        close(4)
      end if

      if(qopt) then
        call openifile(1, 4,chfile)         ! TDS Table file
        do i=1,50
          read(4,'(2i8)',end=411) i1,j1
          if(i1.ne.i) goto 411
          ntab(1,i) = j1
          do j=1,j1
            read(4,'(2f8.0)') table(1,i,j,1),table(1,i,j,2)
          end do
        end do
 411    close(4)


        call openifile(1, 4,chfile)         ! TIN Table file
        do i=1,50
          read(4,'(2i8)',end=412) i1,j1
          if(i1.ne.i) goto 412
          ntab(2,i) = j1
          do j=1,j1
            read(4,'(2f8.0)') table(2,i,j,1),table(2,i,j,2)
          end do
        end do
 412    close(4)
      end if

      call openifile(1, 5, efile)         ! evaporation file

      call openifile(1, 7, prfile)        ! storm predictor file


      read(1,'(a50,i2)') ofile
        open(6,file=ofile,status='UNKNOWN')
        write(6,'(a,a)') ' Main Input File = ',file0

c     read(1,'(a)') ch1
      read(1,'(a50)') divfile
      write(*,'(1x,a)') divfile
      if(divfile.ne.bl50) open(10,file=divfile,status='UNKNOWN')
c
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

c
c     read hydrologic area names from flow data files
c

      read(11,1043) (haname(kk),kk=1,nhat)
 1043 format((8x,25a8))

      read(11,'((8x,25f8.0))') (haarea(kk),kk=1,nhat)

      ! read start date
      read(11,'(i4,2i2)') isyr,ismo,isdy
      backspace(11)
      write(*,'(a,i4.4,2i2.2)') ' Simulation starts at ',isyr,ismo,isdy

      if(nbf.gt.0) then
        read(12,'(8x,30a8)') (haname(kk),kk=nhat+1,nhat1)
        ! position the file at the start of new simulation
        call tsfpos(12,isyr,ismo,isdy,bifile)
      end if

      if(nps.gt.0) then
        read(13,'(8x,50a8)') (haname(kk),kk=nhat1+1,nhat2)
        ! position the file at the start of new simulation
        call tsfpos(13,isyr,ismo,isdy,psfile)
      end if

      ! position the file at the start of simulation
      call tsfpos(5,isyr,ismo,isdy,efile)
      call tsfpos(7,isyr,ismo,isdy,prfile)

      write(6,'(/,a)') ' Source of water ============================='
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
              jerror = jerror + 1
            end if
          end do
        end if
      end do
      write(6,'(a1)') ' '

      ! set up link and node system
      call LNSetup()

      call rdoutp()


      do i=1,nqls          ! number of reservoirs to print monthly percolation 
        do k=1,nres        ! number of reservoirs in the system
          if(lkatrv(k).eq.ipql(i)) then               
            kpql(i) = k    ! found the reservoir for link ID and save reservoir serial number
            goto 58
          end if
        end do
        write(6,'(a,3i5)')' problem in nqls spec',i
 58     continue
      end do
      write(6,'(//,a)') ' Checking pointers'
      do i=1,nqls
         write(6,'(2i6,2x,a8,i6)')
     &                     i,ipql(i),lkname(ipql(i)),kpql(i)
      end do


      write(6,'(//,a)') 'Channel Flow-Width Relaionship'
      do ich=0,MXCH
        mtab = max(ntab(0,ich), ntab(1,ich), ntab(2,ich))
        if(mtab.eq.0) goto 1464
        write(6,'(i10,3(i10,10x))') ich, (ntab(k,ich),k=0,2)
        do i=1,mtab
          write(6,'(6f10.2)') ((table(k,ich,i,j),j=1,2),k=0,2)
        end do
 1464   continue
      end do

      if(ierror) then
        write(*,'(a)') ' IERROR = .TRUE.'    ! ro_v1 3/20/06
        stop                                 ! ro_v1 3/20/06
      end if                                 ! ro_v1 3/20/06

c
c     start routing by calling each link one at a time and storing
c     one time period at a time
c
      t=0
 1042 format(a1)

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



c
c     ============================= Time Loop ==========================
c
 100  t=t+1

      qout = 0.             ! (MXLK,0:MXPS)     flow out of nodes
      qlout = 0.            ! (MXLK,0:MXPS)     flow at end of links
      cqout = 0.            ! (MXLK,4,0:MXPS)   concentration out of nodes
      cqlout = 0.           ! (MXLK,4,0:MXPS)   concentration at end of links
      qloss = 0.            ! (MXLK,0:MXPS)     percolation from links
      qevp = 0.             ! (MXLK)            evaporation from links
c
      read(11,1006,end=999) iyr,imo,idy,(qha(k),k=1,nhat)
      if(imo.eq.ismo .and. idy.eq.1) then     
        write(*,'(i5,2i2)') iyr               
c        ieyr = iyr                           
      end if                                  
c     if(idy.eq.1) write(*,'(i5,2i2)') iyr,imo  !,idy
c     write(*,'(1x,3i4,a)')iyr,imo,idy, '  ====' ! delete
      mday = mndays(iyr,imo)
        ieyr = iyr+1          

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
          ! calculate TDS and TIN concentration, only option 3 is allowed now
          do k=nhat+1,nhat1
            cqha(k,1) = 0.
            cqha(k,2) = 0.
            if(qha(k).gt.0.) then
              ! TDS
              flow = qha(k) * flowconv(k,1)
              cqha(k,1) = rintp(1,ptrc(k,1),flow)
              if(ldebug) write(*,'(a)') ' pass q0a2'

              ! TIN
              if(ldebug) write(*,'(a)') ' pass q0a2'
              flow = qha(k) * flowconv(k,2)
              cqha(k,2) = rintp(2,ptrc(k,2),flow)
              if(ldebug) write(*,'(a)') ' pass q0a3'
            end if
          end do
        end if
      end if
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

      if(ldebug) write(*,'(a)') ' pass q1b'

      read(5,'(i4,2i2,10f8.0)',end=999) jyr,jmo,jdy,
     &                                  (evapor(k),k=1,newts)
      call cmpdate(iyr,imo,idy,  5,jyr,jmo,jdy)


 1006 format(i4,2i2,25f8.0,/, (8x, 25f8.0))


      do 500 j=1,nlinks ! ==============================================

      if(ldebug) write(*,'(a,i5,1x,a8)') '  Pass q2, j = ',j,lkname(j)
      if(convtype(j).eq.5) then
c
c       dummy split link
c       do nothing
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
        rval1 = qout(iun,0) / rval0
        do kps=1,nps
          qout(iun,kps) = qout(iun,kps) * rval1
        end do
      end if

      if(ldebug) write(*,'(a,2i5)') '  Pass q2c',j,convtype(j)

      goto (150,155,200,130,400,140) convtype(j)  ! ====================

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
      width(j) = rintp(0,ich,flow)
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

      stor0 = storl(k,0) + qin_af                        ! save for concentration calculation
      cstor1 = 0.
      cstor2 = 0.
      if(qopt .and. stor0.gt.0.) then
        cstor1 = (storl(k,0) * cstorl(k,1,0)             ! initial concentration
     &            + qin_af * cqout(iun,1,0)) / stor0
        cstor2 = (storl(k,0) * cstorl(k,2,0)
     &            + qin_af * cqout(iun,2,0)) / stor0
      end if


      call rvsim2(k,qin_af)                ! rvsim2 has storm mode operation

      r1 = 0.
      r2 = 0.
      r3 = 0.
      r4 = 0.
      r5 = 0.
      r6 = 0.
      r7 = 0.
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
      end if

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
          xqloss(i) = xqloss(i) + qloss(ipql(i),0) * cf1
          xqevp(i) = xqevp(i) + qevp(ipql(i))
        end do
        if(idy.eq.mday) then
          write(81,'(i4,i2.2,(50f8.2))') iyr,imo,
     &                              (xqloss(i),i=1,nqls)
          write(82,'(i4,i2.2,(50f8.2))') iyr,imo,
     &                              (xqevp(i),i=1,nqls)
          do i=1,nqls
            xqloss(i) = 0.
            xqevp(i) = 0.
          end do
        end if
        inquire(unit=83, opened = okfile)                         !v6a1
        if(okfile) write(83,'(i4,2i2.2,(50f8.2))') iyr,imo,idy,   !v6a1
     &                  (storl(kpql(i),0),i=1,nqls)               !v6a1
      end if
      if(ldebug) write(*,'(a)') '  Pass q6c' !delete ---------------------------

      if(nprch.gt.0) then
        do j=1,nlinks
          ii = linkrch(j)
          if(ii.gt.0 .and. ii.le.30) then
            iun = usn(j)
            do kps=0,nps
              archp(ii,0,kps,iyr,imo) = archp(ii,0,kps,iyr,imo)
     &                        + qloss(j,kps)
              if(qopt) archp(ii,1,kps,iyr,imo) = archp(ii,1,kps,iyr,imo)
     &                        + qloss(j,kps) * cqlout(j,1,kps)
              if(qopt) archp(ii,2,kps,iyr,imo) = archp(ii,2,kps,iyr,imo)
     &                        + qloss(j,kps)
     &                        * (cqlout(j,2,kps) + cqout(iun,2,kps))/2.
            end do
          end if
        end do
      end if
      if(ldebug) write(*,'(a)') '  Pass q6d' !delete ---------------------------

      if(npnd.gt.0) then
        do i=1,npnd
          ii = ipnd(i)
          do kps=0,nps
            aflow(i,0,kps,iyr,imo) = aflow(i,0,kps,iyr,imo)
     &                             + qout(ii,kps)
            if(qopt) aflow(i,1,kps,iyr,imo) = aflow(i,1,kps,iyr,imo)
     &                             + qout(ii,kps) * cqout(ii,1,kps)
            if(qopt) aflow(i,2,kps,iyr,imo) = aflow(i,2,kps,iyr,imo)
     &                              + qout(ii,kps) * cqout(ii,2,kps)
          end do
        end do
      end if !npnd
      if(ldebug) write(*,'(a)') '  Pass 6e' !delete ---------------------------

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
cc    do iyr=1949,2006
c     do iyr=1994,2006
c     do iyr=1994,2009
      do iyr=1920,2011
        do imo=1,12
          do kps=0,nps
            do ii=1,npnd
              aflow(ii,0,kps,iyr,imo) = aflow(ii,0,kps,iyr,imo) * cf1
              aflow(ii,1,kps,iyr,imo) = aflow(ii,1,kps,iyr,imo) * cf2
              aflow(ii,2,kps,iyr,imo) = aflow(ii,2,kps,iyr,imo) * cf2
            end do
          end do
        end do
      end do

c       write(64,'(a)') ' Total monthly flow passing the nodes'
c       write(64,'(a)') ' Units (ac-ft)'

        kps = 0
c        do ify=1921,2011
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
     &           icy,icm,(aflow(ii,0,kps,icy,icm),ii=1,npnd)
          end do
        end do
      end if


      if(nqls.gt.0) close(81)
      if(nqls.gt.0) close(82)

      if(nprch.gt.0) then !---------------------------------------------
        do iyr=isyr,ieyr
          do imo=1,12
            do kps=0,nps
              do ii=1,nprch
                archp(ii,0,kps,iyr,imo) = archp(ii,0,kps,iyr,imo) * cf1
                archp(ii,1,kps,iyr,imo) = archp(ii,1,kps,iyr,imo) * cf2
                archp(ii,2,kps,iyr,imo) = archp(ii,2,kps,iyr,imo) * cf2
              end do
            end do
          end do
        end do 


        if(qopt) write(68,'(a)') ' Total monthly TDS mass percolation'
        if(qopt) write(68,'(a)') ' Units|(tons)'
        if(qopt) write(69,'(a)') ' Total monthly TIN mass percolation'
        if(qopt) write(69,'(a)') ' Units|(tons)'

        kps = 0
        do ify=isyr,ieyr
          do ifm=1,12
            if(ifm.le.6) then
              icy = ify - 1
              icm = ifm+6
            else
              icy = ify
              icm = ifm-6
            end if
            write(67,'(i4,i2,30f8.1)')
     &           icy,icm,(archp(ii,0,kps,icy,icm),ii=1,nprch)
          end do
        end do
      end if                                                      !?

      close(67)
      if(qopt) close(68)
      if(qopt) close(69)

      stop
      end




