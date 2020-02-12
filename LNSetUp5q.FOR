c LNSetup5q.FOR is modified from LNSetUp5.FOR  2/6/2018
c It includes planned diversion from reservoir/recharge basins
c
c Data input for link types
c link type 1
c
c
c
c

      subroutine LNSetup()

      include 'RO_DIM_3.MAX'

c---- general model control properties
      include 'GenParm.VAR'

c---- link/node properties
      include 'LKprpt2.var'

c---- HSA properties
      include 'HAprpt.var'

c---- reservoir properties
      include 'RVprpt.VAR'

c---- diversion properties
      include 'DVprpt.VAR'

c---- various table properties
      include 'TBprpt.VAR'

      character*8  ch8(MXHN), bl8
      real         rtemp(MXHA)


      bl8 = '        '
      velocity = 1.0
c
c     define node properties
c
      write(*,'(/,a)') ' Reading node data'
      j=0
 10   j=j+1
 20   read(3,'(a8,2i8)') ndname(j),nha(j)
      if(ndname(j).eq.bl8) goto 20      
      call ctrim8r(ndname(j))
      if(ldebug) write(*,'(i10,2x,a8)') j,ndname(j)
      if(ndname(j).eq.'END     ') then
        nnodes = j-1
        goto 40
      end if

      if(nha(j).le.0) go to 10
      if(nha(j).gt.MXHN) then
        write(*,'(a,i4,a10,a)') '  Number of loads at node',j,ndname(j),
     &                         ' is greater than MXHN'
        stop
      end if
      read(3,'(10a8)') (ch8(k),k=1,nha(j))
      do k=1,nha(j)
        call ctrim8r(ch8(k))
      end do
      ! read multiplier
      ! this is for small mountain watershed
      read(3,'(10f8.0)') (rtemp(k),k=1,nha(j))

      do k=1,nha(j)
        do kk=1,nhat2
          if(ch8(k).eq.haname(kk)) then
            nodha(j,k) = kk
            rmult(j,k) = rtemp(k)
            go to 30
          end if
        end do
        write(6,'(a,i2,1x,a8,a,i4,1x,a8,a)') ' HA name ',k,ch8(k),
     &           ' for node #',j,ndname(j),' is not available'
        ierror = .true.
 30     continue
      end do
      go to 10

 40   close(3)

c
c     read link data
c
      write(*,'(/,a)') ' Reading link data'

      nres = 0   ! nres = reservoir counter
      ndiv = 0   ! ndiv = diversion counter
      j = 0      ! j    = node counter

 50   j = j+1
 60   read(2,'(3a8,5i8)') lkname(j),ch8(1),ch8(2),convtype(j)
      if(lkname(j).eq.bl8) goto 60
      call ctrim8r(lkname(j))
      write(*,*) j,lkname(j),ch8(1),ch8(2),convtype(j)
      if(lkname(j).eq.'END     ') goto 600

      call ctrim8r(ch8(1))
      call ctrim8r(ch8(2))
      do in=1,nnodes
        if(ch8(1).eq.ndname(in)) then
          usn(j) = in
          goto 70
        end if
      end do
      write(6,'(a,a,a,a,a)')  ' Upstream node   ',ch8(1), ' for link ',
     &          lkname(j),' is not defined.'
      ierror = .true.
 70   do in=1,nnodes
        if(ch8(2).eq.ndname(in)) then
          dsn(j) = in
          goto 80
        end if
      end do
      write(6,'(a,a,a,a,a)')  ' Downstream node ',ch8(2), ' for link ',
     &          lkname(j),' is not defined.'
      ierror = .true.
 80   continue

c
c    convtype = 1 or 2   conveyance only
c             = 3   storage and conveyance
c             = 4   diversion
c             = 5   dummy link
c

      goto (100,200,300,400,500) convtype(j)
        write(6,'(a,i3,a)') ' Link type ',convtype(j),
     &                      ' is not available.  Stop execution.'
        write(6,'(5x,3a8,5i8)') lkname(j),ch8(1),ch8(2),convtype(j)
        ierror = .true.
        goto 500

c     open channel link
 100  read(2,'(a8,10f8.0)') ch8(1),(rtemp(k),k=1,9)
      if(rtemp(8).ne.0) rxn(j) = rtemp(8)
      call ctrim8r(ch8(1))
      if(ch8(1).eq.bl8) then                  ! this is lined channel
        slop(j)    = 0.
        mannn(j)   = 0.
        b(j)       = 0.
        z(j)       = 0.
        length(j)  = 0.
        dpb(j)     = 0.
        dpz(j)     = 0.
      elseif(ch8(1).eq.'MANNING') then        ! this is pervious open channel
        slop(j)    = rtemp(1)                 ! main channel slope
        mannn(j)   = rtemp(2)                 ! manning roughness
        b(j)       = rtemp(3)                 ! trapezoidal bottom width
        z(j)       = rtemp(4)                 ! side slope, assume left and right side slopes are same
        length(j)  = rtemp(5)                 ! channel length
        lkrch(j)   = rtemp(6)                 ! link reach ID to assign parameters from parameter file
        dpb(j)     = rvrch(lkrch(j),1)        ! stream bottom infiltration rate (ft/day)
        dpz(j)     = rvrch(lkrch(j),2)        ! stream side slope infiltration rate - not implemented yet.
        if(dpz(j) .eq. 0.) dpz(j) = dpb(j)    
      else                                    ! this could be flow-width rating option channel
        do ich=1,nch
          if(ch8(1).eq.chname(ich)) then
            convtype(j) = 6
            ptrch(j) = ich
            length(j) = rtemp(5)
            lkrch(j)   = rtemp(6)
            dpb(j)     = rvrch(lkrch(j),1)
            dpz(j)     = rvrch(lkrch(j),2)
            if(dpz(j) .eq. 0.) dpz(j) = dpb(j)
            goto 110
          end if
        end do
        write(*,'(a,1x,a8,a)') ' CH name ',ch8(k),
     &           ' is not available'
        ierror = .true.
        write(*,'(i7,2x,a8,a,a8)') j,lkname(j),' chname ',ch8(k) !!! delete
      end if

 110  rxn(j) = rvrch(lkrch(j),3)
      rxn(j) = rxn(j)*length(j)/(velocity*86400.)

      if(rtemp(9).gt.0. .and. rtemp(9).le.MXRH) linkrch(j) = rtemp(9)    ! this is reach section of the stream, such as RB's SAR Reaches

      ! unit conversion
      dpb(j) = dpb(j)/86400.               ! (ft/day) to (ft/sec)
      dpz(j) = dpz(j)/86400.
      go to 500
 200  read(2,'(8x,10f8.0)') slop(j),mannn(j),b(j)   ! changed 3/22/15
      go to 500

c
c     Link is a reservoir
c
 300  nres=nres+1
      ptrrv(j)=nres
      read(2,'(a40)') rvname(nres)
      write(*,'(16x,i3,1x,a)') nres,rvname(nres)
c      read(2,'(i8)') npts(nres)
      read(2,'(2i8)') npts(nres),kpdiv(j)
	  if(kpdiv(j).gt.npdd) then   
	    write(*,'(a)') ' Problem w/ predefined diversion'
		write(*,'(1x,a8,2i8)') lkname(j),kpdiv(j),npdd
		stop
	  end if

      slimit(nres,1) = 0.                  ! dead storage
      slimit(nres,2) = 0.                  ! storage at spillway starts
      slimit(nres,3) = 0.                  ! maximum storage
      do kk=1,npts(nres)
        read(2,'(10f8.0)') rvelev(nres,kk),rvarea(nres,kk),
     &               rvstor(nres,kk), rvout(nres,kk,1), 
     &               rvout(nres,kk,2),rvpr(nres,kk)
        if(rvout(nres,kk,1).eq.0. .and. rvout(nres,kk,2).eq.0)
     &               slimit(nres,1) = rvstor(nres,kk)
        if(rvout(nres,kk,2).eq.0)
     &               slimit(nres,2) = rvstor(nres,kk)
      end do
      alimit(nres)   = rvarea(nres,npts(nres))
      slimit(nres,3) = rvstor(nres,npts(nres))

      ! check for storage increase
      do kk=1,npts(nres)-1
        if(rvstor(nres,kk).ge.rvstor(nres,kk+1)) then
          write(*,'(a)') ' Problem with reservoir storage'
          write(*,'(a,i3)') ' Reservoir ',nres
          k=nres
          write(*,'(6f8.2)') rvelev(k,kk),rvarea(k,kk),rvstor(k,kk)
          write(*,'(6f8.2)') rvelev(k,kk+1),rvarea(k,kk+1),
     &                       rvstor(k,kk+1)
          stop
        end if
      end do

c     read downstream link for each outflow
      read(2,'(24x,2a8)') rdlink(nres,1), rdlink(nres,2)
      call ctrim8r(rdlink(nres,1))
      call ctrim8r(rdlink(nres,2))
      lkatrv(nres) = j

      read(2,'(10f8.0)') (evwts(nres,kk),kk=1,newts)
c     read(2,'(f8.0,56x,f8.0)') dperc,rxn(j)
      read(2,'(10f8.0)') (rtemp(k),k=1,10)
      rxn(j) = rtemp(9)
      if(rtemp(10).gt.0..and.rtemp(10).le.MXRH) linkrch(j) = rtemp(10)
      krvtype(nres) = rtemp(2)
      ffflow(nres) = rtemp(4)  ! this is the first flush to be wasted on ffmode day

      if(krvtype(nres).lt.1 .or. krvtype(nres).gt.3) then
        write(*,'(a,2i5)') ' Reservoir type for ',nres,krvtype(nres)
        stop
      end if
      slimit(nres,4) = rtemp(3)

c       write(6,'(i4,4f8.1,a)') nres,(rtemp(jjj),jjj=1,4),' rtemp' !delete
c       write(6,'(i4,f8.1,a)') nres,slimit(nres,4), ' slimit(,4)' !delete

      if(rtemp(3).eq.0.) then ! operating storage is not specified
c       write(6,'(i4,f8.1,a)') nres,rtemp(3), ' rtemp(3) ' !delete
c       write(6,'(i4,f8.1,a)') nres,slimit(nres,4), ' slimit(,4)' !delete
        if(krvtype(nres).eq.1)  then
          slimit(nres,4) = slimit(nres,2)
        elseif (krvtype(nres).eq.2)  then
          slimit(nres,4) = slimit(nres,1)
        else
          slimit(nres,4) = (slimit(nres,1) + slimit(nres,2)) / 2.
        end if
      end if
      ! krvtype() = 1 --- conservation basin
      ! krvtype() = 2 --- flood control basin
      ! krvtype() = 3 --- multipurpose basin
      ! for conservation basins, outlet gates are set
      !   for best conservation 5
      ! for flood control basins, outlet gates should be set
      !   to drain reservoir within 24 hours
      ! for multipurpose basins, all gates should be set for conservation
      !   when a storm is coming, all water should be just released
      !   to the first downstream link.

      do kk=1,npts(nres)
        if(rvpr(nres,kk).eq.0.) rvpr(nres,kk) = rtemp(1)
      end do

      ! setup drain flow tables
      k = nres
      call EASSET(k)


      goto 500
c
c     diversion link
c
 400  ndiv = ndiv+1
      ptrdv(j)=ndiv
      read(2,'(a40)') dvname(ndiv)
      read(2,'(i8)') ndpts(ndiv)
      do kk=1,ndpts(ndiv)
        read(2,'(10f8.0)') dflow(ndiv,kk),divout(ndiv,kk)
      end do
c
c     read downstream links
c       linkdv(ndiv,1) is the main link
c       linkdv(ndiv,2) is the diversion link
c
      read(2,'(2a8)') ddlink(ndiv,2),ddlink(ndiv,1)
      call ctrim8r(ddlink(ndiv,1))
      call ctrim8r(ddlink(ndiv,2))
      lkatdv(ndiv) = j
       ! *** diversion can not be made to upstream link
c     goto 500

 500   go to 50
 600   nlinks=j-1
      close(2)

	  write(*,*) ' Now checking diverion links'
c
c     assign link indices to diversion destination links
c
      do k=1,ndiv
        do i2=1,2
          linkdv(k,i2)=0
          do jk=lkatdv(k)+1,nlinks          ! diversion should not be made to a upstream link
           if(ddlink(k,i2).eq.lkname(jk)) then
             linkdv(k,i2)=jk
             go to 610
           end if
          end do
 610      continue
        end do

        if(linkdv(k,1).eq.0 .or. linkdv(k,2).eq.0) then
          write (*,'(a,i3,a)') ' Outflow links for diversion',k,
     &                            ' are not properly specified'
          ierror = .true.
        end if
      end do
	  
	  write(*,*) ' Now checking reservoir down links'
c
c     assign link indices to reservoir outlet links
c
      do k=1,nres
        do i2=1,2
          linkrv(k,i2)=0
          if(rdlink(k,i2).ne.bl8) then
            do jk=lkatrv(k)+1,nlinks
             if(rdlink(k,i2).eq.lkname(jk)) then
               linkrv(k,i2)=jk
               go to 620
             end if
            end do
 620        continue
          end if
        end do

        if( (rdlink(k,1).ne.bl8 .and. linkrv(k,1).eq.0) .or.
     &      (rdlink(k,2).ne.bl8 .and. linkrv(k,2).eq.0) ) then
          write (*,'(a,i3,a)') ' Outflow links for reservoir',k,
     &                         ' is not properly specified'
          ierror = .true.
        end if
      end do

	  write(*,*) ' Now configuring link-node system'
c     ==================================================================
c
c     Configure Link-Node system
c
      nincl(1) = 0
      do j=2,nlinks
        kk=0
        do jj=1,nlinks
          if(usn(j).eq.dsn(JJ)) then
            kk=kk+1
            incoml(usn(j),kk)=jj
          end if
        end do
        nincl(usn(j))=kk
      end do

	  write(*,*) ' Now writing link-node to main output'
	  
c     ==================================================================
c
c     write link-node configuration to main output file
c
      write(6,'(//,a,/)') '   *** Link - Node Configuration ***'
      write(6,'(a43,a21,a)')
     &        'Link Name    Connecting Nodes        Type  ',
     &        'Node Name   NINCL NHA',
     &        '  Incomming Links'
      write(6,'(13x,a26,27x,a)')
     &        'Upstream     Downstream   ',
     &        'HA Column # and Names'

      do j=1,nlinks
	    write(*,*) j
        ! upstream node
        iun = usn(j)
        if(nincl(iun).eq.0) then
          write(6,'(/,43x,i4,1x,a8,2i4,2x,a)')
     &             iun,ndname(iun),nincl(iun),nha(iun),'Headwater'
        else
          write(6,'(43x,i4,1x,a8,2i4,2x,5(i4,1x,a8))')
     &             iun,ndname(iun),nincl(iun),nha(iun),
     &             (incoml(iun,kk),
     &              lkname(incoml(iun,kk)),kk=1,nincl(iun))
        end if
        if(nha(iun).gt.0) then
          write(6,'(66x,10(i4,1x,a8))')
     &             (nodha(iun,kk),haname(nodha(iun,kk)),kk=1,nha(iun))
          write(6,'(66x,10(5x,f8.3))')
     &             (rmult(iun,kk),kk=1,nha(iun))
        end if

        ! write link information
        write(6,'(3(i4,1x,a8),i2)') j,lkname(j),
     &          usn(j),ndname(usn(j)),dsn(j),ndname(dsn(j)),
     &          convtype(j)
        if(convtype(j).eq.3) then
          ires = ptrrv(j)
          write(6,'(5x,a40)') rvname(ires)
          write(6,'(5x,a,2(i4,1x,a8))') 'Downlink',
     &            linkrv(ires,1),lkname(linkrv(ires,1)),
     &            linkrv(ires,2),lkname(linkrv(ires,2))
        elseif(convtype(j).eq.4) then
          idiv = ptrdv(j)
          write(6,'(5x,a40)') dvname(idiv)
          write(6,'(5x,a,2(i4,1x,a8))') 'Downlink',
     &            linkdv(idiv,2),lkname(linkdv(idiv,2)),
     &            linkdv(idiv,1),lkname(linkdv(idiv,1))
        end if
      end do ! nlinks

      ! last node
      idn = dsn(nlinks)
      write(6,'(43x,i4,1x,a8,2i4,2x,5(i4,1x,a8))')
     &         idn,ndname(idn),nincl(idn),nha(idn),
     &         (incoml(idn,kk),lkname(incoml(idn,kk)),kk=1,nincl(idn))
      if(nha(idn).gt.0) then
        write(6,'(66x,10(i4,1x,a8))')
     &           (nodha(idn,kk),haname(nodha(idn,kk)),kk=1,nha(idn))
        write(6,'(66x,10(5x,f8.3))')
     &           (rmult(idn,kk),kk=1,nha(idn))
      end if
      write(6,'(/,a,/)') 'End of Link-Node Configuration'
      write(6,'(a,i5)') 'Number of Links, nlinks = ',nlinks
      write(6,'(a,i5)') 'Number of Nodes, nnodes = ',nnodes


      write(6,'(/,a)') 'Diversions in the system'
      do k=1,ndiv
        write(6,'(2i5,2x,a,2f10.1)') k,lkatdv(k),dvname(k)
      end do
      write(6,'(/,a)') 'Reservoirs in the system'
      do k=1,nres
        write(6,'(2i5,a10,a,5f10.1)') k,lkatrv(k),lkname(lkatrv(k)),
     &               rvname(k),alimit(k),slimit(k,1),slimit(k,2),
     &               slimit(k,3)
      end do

      return
      end
