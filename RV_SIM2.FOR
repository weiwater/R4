

!===========================================================================
      subroutine RVSIM2(k,qin_af)


      include 'ro_dim_3.max'
      include 'genparm.var'
      include 'rvprpt.var'

      cf1 = 86400./43560.

      dpercl = 0.
      devapl = 0.
      dout1l = 0.
      dout2l = 0.


c     adate = date - 1.     !!!!!!!!!!!!!!!!!

      ! call EASOP() to set drain flow rate
      call EASOP(k)

      stor1 = storl(k,0) + qin_af
      if(stor1.le.0.01) then          ! dry basin
        storl(k,0) = stor1
        goto 300
      end if
      devrate = evapl(k)              ! r3 evapl(k) is used to calculate
                                      ! and transfer daily evaporation rate


      ! first flush modification
      if(ffmode) then
        ffbypass = min(qin_af,ffflow(k))
        qin_af = qin_af - ffbypass
      end if

      maxts = 240 ! 1/10 of an hour
      fmin = slimit(k,3)/maxts

c     npd = min(maxts,max(1,int(qin_af/fmin)))
      npd = min(maxts,max(1,int((qin_af+storl(k,0))/fmin)))
      delt = 1./npd                 ! time step for daily simulaiton
      ainflow = qin_af*delt
      stor1 = storl(k,0)


      do ipd=1,npd

        stor2 = stor1 + ainflow/2.         ! save for concentration calculation
        ! get elevation, surface area and outflow rate
        do kk=2,npts(k)-1
          if(stor2.le.rvstor(k,kk)) go to 215
        end do
 215    slopei=(stor2-rvstor(k,kk-1))/(rvstor(k,kk)-rvstor(k,kk-1))
        aarea = rvarea(k,kk-1)+(rvarea(k,kk)-rvarea(k,kk-1))*slopei
        dperc =   rvpr(k,kk-1)+(  rvpr(k,kk)-  rvpr(k,kk-1))*slopei   
        apercl = aarea*dperc*delt
        aevapl = aarea*devrate*delt

        aout1l = rvout(k,kk-1,1)+(rvout(k,kk,1)-rvout(k,kk-1,1))*slopei
        aout2l = rvout(k,kk-1,2)+(rvout(k,kk,2)-rvout(k,kk-1,2))*slopei
        aout1l = aout1l*delt*cf1
        aout2l = aout2l*delt*cf1

        tout = apercl + aevapl + aout1l + aout2l
        stor2 = stor1 + ainflow - tout
        if(stor2.lt.0.) then
          ratio = (stor1 + ainflow) / tout
          tout = stor1 + ainflow
          apercl = apercl * ratio
          aevapl = aevapl * ratio
          aout1l = aout1l * ratio
          aout2l = aout2l * ratio
          stor2 = 0.
        end if

        ! update
        stor1 = stor2
        ! tally for the day
        dpercl = dpercl + apercl
        devapl = devapl + aevapl
        dout1l = dout1l + aout1l
        dout2l = dout2l + aout2l

      end do  ! ipd

      storl(k,0) = stor1
 300  continue

      if(ffmode) dout2l = dout2l + ffbypass
      ! note water quality portion is not updated for this by-pass

      percl(k) = dpercl
      evapl(k) = devapl
      out1l(k) = dout1l / cf1
      out2l(k) = dout2l / cf1

c     write(6,'(8f8.2)')
c    &     date, qin(k), qin_af,dpercl,devapl,dout1l,dout2l,storl(k,0)


      return
      end


!==========================================================================
      ! EASOP() should be calt(k,kk,4)   ! to drain water

      subroutine EASOP(k)

      include 'genparm.var'
      include 'RO_Dim_3.max'
      include 'rvprpt.var'

      ! choose proper drain flow rate for reservoir type and reservoir operation mode

      do kk=1,npts(k)
        if(nrvopmd.ge.2) then 
          rvout(k,kk,1) = rvout(k,kk,4)   ! to drain water
        else
          rvout(k,kk,1) = rvout(k,kk,3)   ! conservation mode
        end if
          
      end do
      return
      end




!==========================================================================
      ! EASSET() needs to be called once for each recharge basins

      subroutine EASSET(k)

      include 'genparm.var'
      include 'RO_Dim_3.max'
      include 'rvprpt.var'

      ! adjust the rating table for operation storage
      ! check operation storage
      if(slimit(k,4).lt.slimit(k,1)) then
        if(ldebug) then
          write(*,'(a,i4,2x,a)') ' Operation storage for',nres,
     &               rvname(k)
          write(*,'(a)') ' will be modfied for dead storage'
        end if
        slimit(k,4) = slimit(k,1)
      elseif(slimit(k,4).gt.slimit(k,2)) then
        if(ldebug) then
          write(*,'(a,i4,2x,a)') ' Operation storage for',nres,
     &               rvname(k)
          write(*,'(a)') ' will be modfied for spillway storage'
        end if
c       slimit(k,4) = slimit(k,2)
c removed for winevile diversion
      end if

cccc  write(*,'(/,i5)') npts(k)
      do kk=1,npts(k)
        if(rvstor(k,kk).le.slimit(k,4)) then
          goto(11,12,13) krvtype(k)
 11       rvout(k,kk,3) = 0.                    ! conservation basin
          rvout(k,kk,4) = 0.
          goto 14
 12       rvout(k,kk,3) = rvout(k,kk,1)         ! flood-control basin
          rvout(k,kk,4) = rvout(k,kk,1)
          goto 14
 13       rvout(k,kk,3) = 0.                    ! multipurpose basin
          rvout(k,kk,4) = rvout(k,kk,1)
        else
          kk1=kk
          goto 15
        end if
 14     continue
      end do
 15   next = 0
      if(rvstor(k,kk1-1).lt.slimit(k,4)) then
        ! need to add an entry
        next = 1
        if(npts(k)+next.gt.MXTB) then
          write(*,'(a,i4,2x,a)') ' Error with res.',nres,rvname(k)
          write(*,'(a)') ' Number of entry > MXTB'
          stop
        end if
      end if
      do kk=npts(k),kk1,-1
        rvelev(k,kk+next) = rvelev(k,kk)
        rvarea(k,kk+next) = rvarea(k,kk)
        rvstor(k,kk+next) = rvstor(k,kk)
        rvout(k,kk+next,1) = rvout(k,kk,1)
        rvout(k,kk+next,2) = rvout(k,kk,2)
        rvout(k,kk+next,3) = rvout(k,kk,1)
        rvout(k,kk+next,4) = rvout(k,kk,1)
        rvpr(k,kk+next) = rvpr(k,kk)
      end do
      if(next.gt.0) then
        r1 = (rvstor(k,kk1+1) - slimit(k,4)) /
     &       (rvstor(k,kk1+1) - rvstor(k,kk1-1))
        r2 = 1. - r1
        rvelev(k,kk1) = rvelev(k,kk1-1)*r1+rvelev(k,kk1+1)*r2
        rvarea(k,kk1) = rvarea(k,kk1-1)*r1+rvarea(k,kk1+1)*r2
        rvstor(k,kk1) = slimit(k,4)
        rvout(k,kk1,1) = rvout(k,kk1-1,1)*r1 +
     &                      rvout(k,kk1+1,1)*r2
        rvout(k,kk1,2) = rvout(k,kk1-1,2)*r1 +
     &                      rvout(k,kk1+1,2)*r2
        rvpr(k,kk1) = rvpr(k,kk1-1)*r1 + rvpr(k,kk1+1)*r2
        goto(21,22,23) krvtype(k)
 21     rvout(k,kk1,3) = 0.
        rvout(k,kk1,4) = 0.
        goto 24
 22     rvout(k,kk1,3) = rvout(k,kk1,1)
        rvout(k,kk1,4) = rvout(k,kk1,1)
        goto 24
 23     rvout(k,kk1,3) = 0.
        rvout(k,kk1,4) = rvout(k,kk1,1)
 24     continue
        npts(k) = npts(k)+next
      end if

      if(ldebug) then
        write(6,'(/,i4,1x,a)') k,rvname(k)
        write(6,'(i8)') npts(k)
        do kk=1,npts(k)
          write(6,'(10f8.1)') rvelev(k,kk),
     &                        rvarea(k,kk),rvstor(k,kk),
     &                        rvout(k,kk,1), rvout(k,kk,2),
     &                        rvpr(k,kk),
     &                        rvout(k,kk,3), rvout(k,kk,4)
        end do
        write(6,'(a8,5f8.1)') 'SLimit  ',(slimit(k,i),i=1,4)
        write(6,'(a,i2)') '   Reservoir type=',krvtype(k)
      end if

      return
      end




!==========================================================================

      subroutine rvopmd(iyr,imn,idy,iu)
      ! this subroutine read storm prediction and set mode of reservoir operation
      ! iyr,imn,idy defines start date of the simulation
      ! nyr         number of years to simulate
      ! iu,ic       rainfall data input file and column number to use
      ! ju          temporary file unit
      ! iu and ju are open already
      ! nrvopmd reservoir operation mode for the day
      !     0   conservation mode
      !     1   pre-storm mode 1 - no supplemental water is accepted to multi-purpose basins
      !     2   pre-storm mode 2 - drain all multi-purpose basins
      !     3   storm mode

      include 'RO_Dim_3.max'
      include 'genparm.var'
      include 'rvprpt.var'

      read(iu,'(i4,2i2,7f8.2)') jyr,jmn,jdy,(rainfc(i),i=1,7)
      call cmpdate(iyr,imn,idy, iu,jyr,jmn,jdy)
      ndaff = ndaff+1
      
      ! first flush operation
      ffmode = .false.
      if(rainfc(1).ge.rainff) then
        if(ndaff.gt.nffday) ffmode = .true.
        ndaff = 0
      end if

      ! reservoir operation mode depending storm forecast
      nrvopmd = 0
      if(rainfc(1).ge.rainss) then
        nrvopmd = 3
      else
        do i=2,3
          if(rainfc(i).ge.rainss) then
            nrvopmd = 2
            goto 10
          end if
        end do
        do i=4,7
          if(rainfc(i).ge.rainss) then
            nrvopmd = 1
            goto 10
          end if
        end do
      end if
 10   continue
      return
      end



!==========================================================================

      subroutine gensf(iyr,imn,idy,nyr,iu,ic,ju)
      ! this subroutine generates a temporary rainfall forecast file
      ! iyr,imn,idy defines start date of the simulation
      ! nyr         number of years to simulate
      ! iu,ic       rainfall data input file and column number to use
      ! ju          temporary file unit
      ! iu and ju are open already

      real rval1(10), rval2(7)

      lyr = iyr+nyr
      lmn = imn
      ldy = idy

      ! locate the file pointer at the start
 1    read(iu,'(i4,2i2)',err=99,end=99) jyr,jmn,jdy
      if(jyr.eq.iyr .and. jmn.eq.imn .and. jdy.eq.idy) then
        backspace(iu)
        goto 2
      end if
      goto 1
 2    do i=1,7
        read(iu,'(i4,2i2,10f8.0)') kyr,kmn,kdy,(rval1(j),j=1,ic)
        rval2(i) = rval1(ic)
      end do
      write(ju,'(i4,2i2,7f8.2)') iyr,imn,idy,(rval2(i),i=1,7)

      do while (.true.)
        do i=1,6
          rval2(i) = rval2(i+1)
        end do
        call adddays(iyr,imn,idy,1, jyr,jmn,jdy)
        iyr=jyr
        imn=jmn
        idy=jdy
        if(iyr.eq.lyr .and. imn.eq.lmn .and. idy.eq.ldy) goto 50

        call adddays(iyr,imn,idy,6, jyr,jmn,jdy)
        read(iu,'(i4,2i2,10f8.0)') kyr,kmn,kdy,(rval1(j),j=1,ic)
        if(kyr.ne.jyr .or. kmn.ne.jmn .or. kdy.ne.jdy) then
          write(*,'(a)') ' Error in rain forecsat data file'
          write(*,'(2(i6,2i2))') jyr,jmn,jdy,kyr,kmn,kdy
          stop
        end if
        rval2(7) = rval1(ic)
        write(ju,'(i4,2i2,7f8.2)') iyr,imn,idy,(rval2(i),i=1,7)
      end do

 50   return

      return
 99   write(*,'(a)') ' Problem with rain forecast data file'
      stop
      end


!==========================================================================
