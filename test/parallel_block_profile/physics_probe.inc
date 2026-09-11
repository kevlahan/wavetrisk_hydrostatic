    subroutine physics_probe(tag)
      ! Diagnostic copy only. Observe existing packed arrays, never repack or
      ! call a numerical producer. Raw solution is still the pre-physics input.
      character(*), intent(in) :: tag
      character(256) :: filename
      character(32) :: setting
      integer :: u,status,selected,l
      call get_environment_variable('WAVETRISK_RESTART_PROBE',setting,status=status)
      if (status/=0.or.trim(setting)/='1') return
      call get_environment_variable('WAVETRISK_PROBE_DOMAIN',setting,status=status)
      selected=-1
      if (status==0.and.len_trim(setting)>0) read(setting,*) selected
      call get_environment_variable('WAVETRISK_PHYSICS_PROBE_ALL',setting,status=status)
      if (status==0.and.trim(setting)=='1') selected=-1
      if (selected>=0.and.glo_id(rank+1,d)/=selected) return
      write(filename,'(a,a,a,i0,a,i0,a)') 'physics-',tag,'-step-',istep_cumul,'-rank-',rank,'.bin'
      open(newunit=u,file=trim(filename),status='unknown',position='append', &
           access='stream',form='unformatted',action='write')
      write(u) 181402,istep_cumul,glo_id(rank+1,d),id,p_null,i,j,is_pole,mask,zlevels,Nsoil
      write(u) time,latitude,longitude,rho_dz
      do l=1,zlevels
         write(u) sol(S_VELO,l)%data(d)%elts(id_edge(id))
         write(u) sol(S_MASS,l)%data(d)%elts(id_i),sol_mean(S_MASS,l)%data(d)%elts(id_i)
         write(u) sol(S_TEMP,l)%data(d)%elts(id_i),sol_mean(S_TEMP,l)%data(d)%elts(id_i)
      end do
      write(u) phys_dt,nth_day,day_fraction,phys_Phisurf,phys_Play,phys_Pint,phys_Phi,phys_Umag, &
           phys_U,phys_V,phys_W,phys_Theta,Tsoil
      close(u)
    end subroutine physics_probe
