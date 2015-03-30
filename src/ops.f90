module ops_mod
  use domain_mod
  use arch_mod
  use viscous_mod
  implicit none

contains
  subroutine init_ops_mod()
      logical :: initialized = .False.
      if (initialized) return ! initialize only once
      call init_domain_mod()
      initialized = .True.
  end subroutine init_ops_mod

  subroutine post_step1(dom, p, c, offs, dims)
      type(Domain) dom
      integer p
      integer c
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,9) :: dims
      integer id, idS, idW, idSW, idN, idE, idNE
      real(8) pv_SW, pv_W, pv_S, pv_LORT, pv_UPLT, pv_SW_LORT, pv_SW_UPLT, pv
      real(8) phi(0:N_BDRY), full_depth(0:N_BDRY)
      phi(0:NORTHEAST) = 1.0

      if (c .eq. IJMINUS) then
          idSW = idx(-1, -1, offs, dims)
          idW = idx(-1, 0, offs, dims)
          idS = idx(0, -1, offs, dims)
          id = idx(0, 0, offs, dims)
          idN = idx(0, 1, offs, dims)
          idE = idx(1, 0, offs, dims)
          if (penalize) phi(0:WEST) = phi(0:WEST) + alpha_m1*penal%data(dom%id+1)%elts((/id,idN,idE,idS,idW/)+1)
          full_depth(0:WEST) = height((/id,idN,idE,idS,idW/)+1) + &
                             dom%topo%elts((/id,idN,idE,idS,idW/)+1) * phi(0:WEST)
          dom%vort%elts(TRIAG*idSW+LORT+1) = &
                      (velo(EDGE*idW+RT+1)*dom%len%elts(EDGE*idW+RT+1) &
                      -velo(EDGE*idSW+1)*dom%len%elts(EDGE*idSW+1) &
                      -velo(EDGE*idS+UP+1)*dom%len%elts(EDGE*idS+UP+1))
          pv_SW = (dom%corolis%elts(TRIAG*idSW+1) + dom%vort%elts(TRIAG*idSW+1))/ &
                       (full_depth(WEST)*dom%areas%elts(idW+1)%part(6) &
                      + full_depth(0)*sum(dom%areas%elts(id+1)%part(4:5)) &
                      + full_depth(SOUTH)*dom%areas%elts(idS+1)%part(3))
          pv_W = (dom%corolis%elts(TRIAG*idW+LORT+1) + dom%vort%elts(TRIAG*idW+LORT+1)*dom%triarea%elts(TRIAG*idW+LORT+1))/ &
                       (full_depth(WEST)*dom%areas%elts(idW+1)%part(1) &
                      + full_depth(0)*dom%areas%elts(id+1)%part(3) &
                      + full_depth(NORTH)*dom%areas%elts(idN+1)%part(5))
          pv_S = (dom%corolis%elts(TRIAG*idS+UPLT+1) + dom%vort%elts(TRIAG*idS+UPLT+1)*dom%triarea%elts(TRIAG*idS+UPLT+1))/ &
                       (full_depth(SOUTH)*dom%areas%elts(idS+1)%part(2) &
                      + full_depth(EAST)*dom%areas%elts(idE+1)%part(4) &
                      + full_depth(0)*dom%areas%elts(id+1)%part(6))
          dom%vort%elts(TRIAG*idSW+LORT+1) = dom%vort%elts(LORT+TRIAG*idSW+1)/dom%triarea%elts(LORT+TRIAG*idSW+1)
          dom%vort%elts(TRIAG*idSW+UPLT+1) = dom%vort%elts(LORT+TRIAG*idSW+1)
          dom%qe%elts(EDGE*idW+RT+1) = 0.5_8*(pv_W + pv_SW)
          dom%qe%elts(EDGE*idS+UP+1) = 0.5_8*(pv_S + pv_SW)
      end if
      if (c .eq. IPLUSJMINUS) then
          id   = idx(PATCH_SIZE,    0, offs, dims)
          idSW = idx(PATCH_SIZE-1, -1, offs, dims)
          idS  = idx(PATCH_SIZE,   -1, offs, dims)
          idW  = idx(PATCH_SIZE-1,  0, offs, dims)
          idE  = idx(PATCH_SIZE+1,  0, offs, dims)
          idNE = idx(PATCH_SIZE+1,  1, offs, dims)
          if (penalize) phi(0:NORTHEAST) = phi(0:NORTHEAST) + alpha_m1*penal%data(dom%id+1)%elts((/id,id,idE,idS,idW,idNE/)+1)
          full_depth(0:NORTHEAST) = height((/id,id,idE,idS,idW,idNE/)+1) + &
                             dom%topo%elts((/id,id,idE,idS,idW,idNE/)+1) * phi(0:NORTHEAST)
          dom%vort%elts(LORT+TRIAG*idSW+1) = - &
                  ((velo(EDGE*idSW+RT+1)*dom%len%elts(EDGE*idSW+RT+1) + &
                  velo(DG+EDGE*idSW+1)*dom%len%elts(DG+EDGE*idSW+1)) - &
                  velo(EDGE*id+RT+1)*dom%len%elts(EDGE*id+RT+1))
          pv_SW_LORT = (dom%corolis%elts(LORT+TRIAG*idSW+1) + &
                 dom%vort%elts(LORT+TRIAG*idSW+1))/( &
                 full_depth(SOUTHWEST)*dom%areas%elts(idSW+1)%part(1) + &
                 full_depth(SOUTH)*dom%areas%elts(idS +1)%part(3) + &
                 full_depth(0)*sum(dom%areas%elts(id+1)%part(5:6)))
          pv_LORT = (dom%corolis%elts(TRIAG*id+LORT+1) + dom%vort%elts(TRIAG*id+LORT+1)*dom%triarea%elts(TRIAG*id+LORT+1))/ &
                       (full_depth(0)*dom%areas%elts(id  +1)%part(1) &
                      + full_depth(EAST)*dom%areas%elts(idE +1)%part(3) &
                      + full_depth(NORTHEAST)*dom%areas%elts(idNE+1)%part(5))
          pv_SW_UPLT = (dom%corolis%elts(TRIAG*idSW+UPLT+1) + &
                        dom%vort%elts(TRIAG*idSW+UPLT+1)*dom%triarea%elts(TRIAG*idSW+UPLT+1))/ &
                       (full_depth(SOUTHWEST)*dom%areas%elts(idSW+1)%part(2) &
                      + full_depth(0)*dom%areas%elts(id  +1)%part(4) &
                      + full_depth(WEST)*dom%areas%elts(idW +1)%part(6))
          dom%vort%elts(LORT+TRIAG*idSW+1) = dom%vort%elts(LORT+TRIAG*idSW+1)/dom%triarea%elts(LORT+TRIAG*idSW+1)
          dom%vort%elts(TRIAG*idS+UPLT+1) = dom%vort%elts(LORT+TRIAG*idSW+1)
          pv_S = pv_SW_LORT
          dom%qe%elts(EDGE*id  +RT+1) = 0.5_8*(pv_S + pv_LORT)
          dom%qe%elts(EDGE*idSW+DG+1) = 0.5_8*(pv_SW_LORT + pv_SW_UPLT)
      end if
      if (c .eq. IMINUSJPLUS) then
          id   = idx(0,  PATCH_SIZE,   offs, dims)
          idSW = idx(-1, PATCH_SIZE-1, offs, dims)
          idW  = idx(-1, PATCH_SIZE,   offs, dims)
          idS  = idx(0,  PATCH_SIZE-1, offs, dims)
          idN  = idx(0,  PATCH_SIZE+1, offs, dims)
          idNE = idx(1,  PATCH_SIZE+1, offs, dims)
          if (penalize) phi(0:NORTHEAST) = phi(0:NORTHEAST) + alpha_m1*penal%data(dom%id+1)%elts((/id,idN,id,idS,idW,idNE/)+1)
          full_depth(0:NORTHEAST) = height((/id,idN,id,idS,idW,idNE/)+1) + &
                             dom%topo%elts((/id,idN,id,idS,idW,idNE/)+1) * phi(0:NORTHEAST)
          dom%vort%elts(TRIAG*idSW+UPLT+1) = &
                  - velo(EDGE*id+UP+1)*dom%len%elts(EDGE*id+UP+1) &
                  + velo(DG+EDGE*idSW+1)*dom%len%elts(DG+EDGE*idSW+1) &
                  + velo(EDGE*idSW+UP+1)*dom%len%elts(EDGE*idSW+UP+1)
          pv_SW_UPLT = (dom%corolis%elts(TRIAG*idSW+UPLT+1) &
                  + dom%vort%elts(TRIAG*idSW+UPLT+1)) &
                  /(full_depth(SOUTHWEST)*dom%areas%elts(idSW  +1)%part(2) &
                  + full_depth(0)*sum(dom%areas%elts(id  +1)%part(3:4)) &
                  + full_depth(WEST)*dom%areas%elts(idW +1)%part(6))
          pv_UPLT = (dom%corolis%elts(TRIAG*id+UPLT+1) + dom%vort%elts(TRIAG*id+UPLT+1)*dom%triarea%elts(TRIAG*id+UPLT+1))/ &
                       (full_depth(0)*dom%areas%elts(id+1)%part(2) &
                      + full_depth(NORTHEAST)*dom%areas%elts(idNE+1)%part(4) &
                      + full_depth(NORTH)*dom%areas%elts(idN+1)%part(6))
          pv_SW_LORT = (dom%corolis%elts(TRIAG*idSW+LORT+1) + &
                        dom%vort%elts(TRIAG*idSW+LORT+1)*dom%triarea%elts(TRIAG*idSW+LORT+1))/ &
                       (full_depth(SOUTHWEST)*dom%areas%elts(idSW+1)%part(1) &
                      + full_depth(SOUTH)*dom%areas%elts(idS+1)%part(3) &
                      + full_depth(0)*dom%areas%elts(id+1)%part(5))
          dom%vort%elts(TRIAG*idSW+UPLT+1) = dom%vort%elts(TRIAG*idSW+UPLT+1)/dom%triarea%elts(TRIAG*idSW+UPLT+1)  
          dom%vort%elts(LORT+TRIAG*idW+1) = dom%vort%elts(TRIAG*idSW+UPLT+1)
          pv_W = pv_SW_UPLT
          dom%qe%elts(EDGE*id  +UP+1) = 0.5_8*(pv_W + pv_UPLT)
          dom%qe%elts(EDGE*idSW+DG+1) = 0.5_8*(pv_SW_LORT + pv_SW_UPLT)
      end if
      if (c .eq. IJPLUS) then
          id  = idx(PATCH_SIZE, PATCH_SIZE, offs, dims)
          idN = idx(PATCH_SIZE, PATCH_SIZE+1, offs, dims)
          idE = idx(PATCH_SIZE+1, PATCH_SIZE, offs, dims)
          idS = idx(PATCH_SIZE, PATCH_SIZE-1, offs, dims)
          idW = idx(PATCH_SIZE-1, PATCH_SIZE, offs, dims)
          if (penalize) phi(0:WEST) = phi(0:WEST) + alpha_m1*penal%data(dom%id+1)%elts((/id,idN,idE,idS,idW/)+1)
          full_depth(0:WEST) = height((/id,idN,idE,idS,idW/)+1) + &
                             dom%topo%elts((/id,idN,idE,idS,idW/)+1) * phi(0:WEST)
          dom%vort%elts(LORT+TRIAG*id+1) = - &
                  (velo(EDGE*id +RT+1)*dom%len%elts(EDGE*id+RT+1) - &
                   velo(EDGE*idN+RT+1)*dom%len%elts(EDGE*id+DG+1) - &
                   velo(EDGE*id +UP+1)*dom%len%elts(EDGE*id+UP+1))
          pv = (dom%corolis%elts(TRIAG*id+1) + dom%vort%elts(LORT+TRIAG*id+1))/ &          
                       (full_depth(EAST)*dom%areas%elts(idE+1)%part(3) + &
                        full_depth(0)*sum(dom%areas%elts(id+1)%part(1:2)) + &
                        full_depth(NORTH)*dom%areas%elts(idN+1)%part(6))
          pv_W = (dom%corolis%elts(TRIAG*idW+LORT+1) + dom%vort%elts(TRIAG*idW+LORT+1)*dom%triarea%elts(TRIAG*idW+LORT+1))/ &
                       (full_depth(WEST)*dom%areas%elts(idW+1)%part(1) &
                      + full_depth(0)*dom%areas%elts(id+1)%part(3) &
                      + full_depth(NORTH)*dom%areas%elts(idN+1)%part(5))
          pv_S = (dom%corolis%elts(TRIAG*idS+UPLT+1) + dom%vort%elts(TRIAG*idS+UPLT+1)*dom%triarea%elts(TRIAG*idS+UPLT+1))/ &
                       (full_depth(SOUTH)*dom%areas%elts(idS+1)%part(2) &
                      + full_depth(EAST)*dom%areas%elts(idE+1)%part(4) &
                      + full_depth(0)*dom%areas%elts(id+1)%part(6))
          dom%vort%elts(LORT+TRIAG*id+1) = dom%vort%elts(LORT+TRIAG*id+1)/dom%triarea%elts(LORT+TRIAG*id+1)
          dom%vort%elts(TRIAG*id+UPLT+1) = dom%vort%elts(LORT+TRIAG*id+1)
          dom%qe%elts(EDGE*id+RT+1) = 0.5_8*(pv + pv_S)
          dom%qe%elts(EDGE*id+UP+1) = 0.5_8*(pv + pv_W)
      end if
  end subroutine

  subroutine step1(dom, p)
      type(Domain) dom
      integer p
      integer i, j, id, offs(0:N_BDRY), dims(2,N_BDRY)
      integer n, e, s, w, ne, sw
      real(8) u_prim_up, u_dual_up, u_prim_dg, u_dual_dg, u_prim_rt, u_dual_rt
      real(8) u_prim_dn, u_dual_dn, u_prim_sw, u_dual_sw, u_prim_lt, u_dual_lt
      real(8) pv_LORT, pv_UPLT, pv_S, pv_W, vort_W, vort_S, vort_LORT, vort_UPLT
      logical S_bdry, W_bdry
      real(8) phi(0:N_BDRY), full_depth(0:N_BDRY)
      
      call comp_offs3(dom, p, offs, dims)
      S_bdry = (dom%patch%elts(p+1)%neigh(SOUTH) .lt. 0)
      if (S_bdry) S_bdry = (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(SOUTH)+1)%side .gt. 0)
      W_bdry = (dom%patch%elts(p+1)%neigh(WEST) .lt. 0)
      if (W_bdry) W_bdry = (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(WEST)+1)%side .gt. 0)

      id = offs(0)
                                     n = PATCH_SIZE  ;  ne = n+1
       w = offs(WEST)                                ;   e =  +1
      sw = offs(SOUTHWEST)        ;  s = offs(SOUTH) 
      call comput()
      if (W_bdry .or. S_bdry) call comp_ijmin()

      w = -1; sw = s+w
      do id = offs(0)+1, offs(0)+LAST-1
          call comput();
          if (S_bdry) call comp_ijmin()
      end do

      e = offs(EAST); ne = dims(1,EAST) + e
      id = offs(0)+LAST
      call comput()
      if (S_bdry) call comp_ijmin()

      s = -PATCH_SIZE
      do j = 2, PATCH_SIZE-1
          id = offs(0) + PATCH_SIZE*(j-1)
          e = +1; ne = n+e
          w = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*(j-1) ! correct for dimension smaller than patch if boundary
          sw = w-dims(1,WEST)
          call comput()
          if (W_bdry) call comp_ijmin()

          w = -1; sw = s+w; ne = n+e
          do id = offs(0)+PATCH_SIZE*(j-1)+1, offs(0)+PATCH_SIZE*(j-1)+LAST-1
              call comput()
          end do
  
          id = offs(0)+PATCH_SIZE*j-1
          e = offs(EAST) + (dims(1,EAST)-PATCH_SIZE)*(j-1); ne = e+dims(1,EAST)
          call comput()
      end do

      id = offs(0)+PATCH_SIZE*LAST
      n = offs(NORTH)
      w = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*LAST
      sw = w-dims(1,WEST)
      e = +1; ne = n+e
      call comput()
      if (W_bdry) call comp_ijmin()
  
      w = -1; sw = s+w
      do id = offs(0)+PATCH_SIZE*LAST+1, offs(0)+PATCH_SIZE*PATCH_SIZE-2
          call comput()
      end do
  
      id = offs(0)+PATCH_SIZE*PATCH_SIZE-1
      ne = offs(NORTHEAST)
      e = offs(EAST) + (dims(1,EAST)-PATCH_SIZE)*LAST
      call comput()

      if (dom%patch%elts(p+1)%neigh(NORTH) .lt. 0) then ! neighbour is boundary
        if (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(NORTH)+1)%side .gt. 0) then ! domain bdry
          id = offs(0)+PATCH_SIZE*LAST + offs(NORTH) ! id + n
           s =                                   - offs(NORTH) ! relative to current id
           w = offs(NORTHWEST)                   - offs(NORTH) 
          sw = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*LAST - offs(NORTH)
           e = 1
           n = dims(1,NORTH)
           ne = n+e
          call comput()
          if (W_bdry) call comp_ijmin()

           w = -1
          sw = -offs(NORTH) + w
          do id = offs(0)+PATCH_SIZE*LAST+1+offs(NORTH), offs(0)+PATCH_SIZE*PATCH_SIZE+offs(NORTH)-2
              call comput()
          end do

           e = offs(NORTHEAST)              - offs(NORTH) 
          ne = e+dims(1,NORTHEAST)
          id = offs(0)+PATCH_SIZE*PATCH_SIZE+offs(NORTH)-1
          call comput()

          ! NORTHEAST point
          id = offs(0)+PATCH_SIZE*PATCH_SIZE+offs(NORTHEAST)-1
           s = -offs(NORTHEAST) + offs(EAST) + (dims(1,EAST)-PATCH_SIZE)*LAST
           sw= -offs(NORTHEAST)
           w = -offs(NORTHEAST) + offs(NORTH)
           n = dims(1,NORTHEAST)
           e = 1
          ne = n+e
          call comput()
        end if
      end if

      if (dom%patch%elts(p+1)%neigh(EAST) .lt. 0) then ! neighbour is boundary
          if (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(EAST)+1)%side .le. 0) return ! not domain bdry
          id = offs(0) + LAST + offs(EAST) ! id + n
          w  =                - offs(EAST) ! relative to current id
          s  = offs(SOUTHEAST)- offs(EAST) 
          sw = offs(SOUTH)    - offs(EAST)
          n  = dims(1,EAST)
          e  = 1
          ne = n+e
          call comput()
          if (S_bdry) call comp_ijmin()
        
          s = -dims(1,EAST)
          do j = 1, LAST-1
              id = offs(0) + LAST + offs(EAST) + j*dims(1,EAST)
              w  =                - offs(EAST) + j*(PATCH_SIZE-dims(1,EAST))
              sw = w-PATCH_SIZE
              call comput()
          end do

          id = offs(0) + LAST + offs(EAST) + LAST*dims(1,EAST)
          n  = offs(NORTHEAST) + LAST*(PATCH_SIZE-dims(1,EAST)) - offs(EAST)
          w  =                - offs(EAST) + LAST*(PATCH_SIZE-dims(1,EAST))
          sw = w-PATCH_SIZE
          ne = n+e
          call comput()
      end if
    contains
      subroutine comp_ijmin()
          real(8) vort_SW
          phi(SOUTHWEST) = 1.0
          if (penalize) phi(SOUTHWEST) = phi(SOUTHWEST) + alpha_m1*penal%data(dom%id+1)%elts(id+sw+1)
          full_depth(SOUTHWEST) = height(id+sw+1) + dom%topo%elts(id+sw+1) * phi(SOUTHWEST)
          vort_SW = - (velo(EDGE*(id+sw)+RT+1)*dom%len%elts(EDGE*(id+sw)+RT+1) + u_prim_sw + u_prim_dn)
          pv_LORT = (dom%corolis%elts(TRIAG*(id+sw)+LORT+1) + vort_SW)/( &
              full_depth(SOUTHWEST)*dom%areas%elts(id+sw+1)%part(1) + &
              full_depth(SOUTH)*dom%areas%elts(id+s +1)%part(3) + &
              full_depth(0)*dom%areas%elts(id   +1)%part(5))
          vort_SW = u_prim_lt + u_prim_sw + velo(EDGE*(id+sw)+UP+1)*dom%len%elts(EDGE*(id+sw)+UP+1) 
          pv_UPLT = (dom%corolis%elts(TRIAG*(id+sw)+UPLT+1) + vort_SW)/( &
              full_depth(SOUTHWEST)*dom%areas%elts(id+sw+1)%part(2) + &
              full_depth(0)*dom%areas%elts(id   +1)%part(4) + &
              full_depth(WEST)*dom%areas%elts(id+w +1)%part(6))
          dom%vort%elts(TRIAG*(id+w)+LORT+1) = vort_W/dom%triarea%elts(TRIAG*(id+w)+LORT+1) 
          dom%vort%elts(TRIAG*(id+s)+UPLT+1) = vort_S/dom%triarea%elts(TRIAG*(id+s)+UPLT+1) 
          dom%qe%elts(EDGE*(id+ w)+RT+1) = 0.5_8*(pv_W   +pv_UPLT)
          dom%qe%elts(EDGE*(id+sw)+DG+1) = 0.5_8*(pv_LORT+pv_UPLT)
          dom%qe%elts(EDGE*(id+s )+UP+1) = 0.5_8*(pv_LORT+pv_S)
          tflux(EDGE*(id+s )+UP+1) = u_dual_dn*(full_depth(SOUTH) + full_depth(0))*0.5_8
          tflux(EDGE*(id+sw)+DG+1) = u_dual_sw*(full_depth(0) + full_depth(SOUTHWEST))*0.5_8
          tflux(EDGE*(id+ w)+RT+1) = u_dual_lt*(full_depth(WEST) + full_depth(0))*0.5_8
      end subroutine
      subroutine comput()
          u_prim_up = velo(EDGE*id+UP+1)*dom%len%elts(EDGE*id+UP+1)
          u_dual_up = velo(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)
          u_prim_dg = velo(EDGE*id+DG+1)*dom%len%elts(EDGE*id+DG+1)
          u_dual_dg = velo(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1)
          u_prim_rt = velo(EDGE*id+RT+1)*dom%len%elts(EDGE*id+RT+1)
          u_dual_rt = velo(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)
          phi(0:NORTHEAST) = 1.0
          if (penalize) phi(0:NORTHEAST) = phi(0:NORTHEAST) + &
                  alpha_m1*penal%data(dom%id+1)%elts(id+(/0,n,e,s,w,ne/)+1)
          full_depth(0:NORTHEAST) = height(id+(/0,n,e,s,w,ne/)+1) + &
                             dom%topo%elts(id+(/0,n,e,s,w,ne/)+1) * phi(0:NORTHEAST)
          tflux(EDGE*id+UP+1) = u_dual_up*(full_depth(0) + full_depth(NORTH))*0.5_8
          tflux(DG+EDGE*id+1) = u_dual_dg*(full_depth(NORTHEAST) + full_depth(0))*0.5_8
          tflux(EDGE*id+RT+1) = u_dual_rt*(full_depth(0) + full_depth(EAST))*0.5_8
          u_prim_dn = velo(EDGE*(id+s)+UP+1)*dom%len%elts(EDGE*(id+s)+UP+1)
          u_dual_dn = velo(EDGE*(id+s)+UP+1)*dom%pedlen%elts(EDGE*(id+s)+UP+1)
          u_prim_sw = velo(EDGE*(id+sw)+DG+1)*dom%len%elts(EDGE*(id+sw)+DG+1)
          u_dual_sw = velo(EDGE*(id+sw)+DG+1)*dom%pedlen%elts(EDGE*(id+sw)+DG+1)
          u_prim_lt = velo(EDGE*(id+w)+RT+1)*dom%len%elts(EDGE*(id+w)+RT+1)
          u_dual_lt = velo(EDGE*(id+W)+RT+1)*dom%pedlen%elts(EDGE*(id+W)+RT+1)
          dom%kin_energy%elts(id+1) = &
              (u_prim_up*u_dual_up + u_prim_dg*u_dual_dg + u_prim_rt*u_dual_rt + &
               u_prim_dn*u_dual_dn + u_prim_sw*u_dual_sw + u_prim_lt*u_dual_lt &
              )* (1.0_8/4.0_8)*dom%areas%elts(id+1)%hex_inv
          if (phi(0) .ne. 0) phi(0) = 1.0_8/phi(0)
          dom%bernoulli%elts(id+1) = grav_accel*height(id+1)*phi(0) + dom%kin_energy%elts(id+1)
          if (viscosity .ne. 0) dom%divu%elts(id+1) = dom%areas%elts(id+1)%hex_inv * &
                  (u_dual_up - u_dual_dg + u_dual_rt - u_dual_dn + u_dual_sw - u_dual_lt)

          dom%vort%elts(LORT+TRIAG*id+1) = - ( &
              u_prim_rt + u_prim_dg + velo(EDGE*(id+E)+UP+1)*dom%len%elts(EDGE*(id+E)+UP+1))
          dom%vort%elts(TRIAG*id+UPLT+1) = &
              velo(EDGE*(id+N)+RT+1)*dom%len%elts(EDGE*(id+N)+RT+1) + u_prim_dg + u_prim_up 
          vort_W = - ( &
              u_prim_lt + velo(EDGE*(id+W)+DG+1)*dom%len%elts(EDGE*(id+W)+DG+1) + u_prim_up)
          vort_S = &
              u_prim_rt + velo(EDGE*(id+S)+DG+1)*dom%len%elts(EDGE*(id+S)+DG+1) + u_prim_dn
          pv_LORT = (dom%corolis%elts(TRIAG*id+LORT+1) + dom%vort%elts(TRIAG*id+LORT+1))/( &
              full_depth(0)*dom%areas%elts(id   +1)%part(1) + &
              full_depth(EAST)*dom%areas%elts(id+E +1)%part(3) + &
              full_depth(NORTHEAST)*dom%areas%elts(id+NE+1)%part(5))
          pv_UPLT = (dom%corolis%elts(TRIAG*id+UPLT+1) + dom%vort%elts(TRIAG*id+UPLT+1))/( &
              full_depth(0)*dom%areas%elts(id   +1)%part(2) + &
              full_depth(NORTHEAST)*dom%areas%elts(id+NE+1)%part(4) + &
              full_depth(NORTH)*dom%areas%elts(id+N +1)%part(6))
          pv_W = (dom%corolis%elts(TRIAG*(id+W)+LORT+1) + vort_W)/( &
              full_depth(WEST)*dom%areas%elts(id+W+1)%part(1) + &
              full_depth(0)*dom%areas%elts(id  +1)%part(3) + &
              full_depth(NORTH)*dom%areas%elts(id+N+1)%part(5))
          pv_S = (dom%corolis%elts(TRIAG*(id+S)+UPLT+1) + vort_S)/( &
              full_depth(SOUTH)*dom%areas%elts(id+S+1)%part(2) + &
              full_depth(EAST)*dom%areas%elts(id+E+1)%part(4) + &
              full_depth(0)*dom%areas%elts(id  +1)%part(6))
          dom%vort%elts(LORT+TRIAG*id+1) = dom%vort%elts(LORT+TRIAG*id+1)/dom%triarea%elts(LORT+TRIAG*id+1) 
          dom%vort%elts(TRIAG*id+UPLT+1) = dom%vort%elts(TRIAG*id+UPLT+1)/dom%triarea%elts(TRIAG*id+UPLT+1) 
          dom%qe%elts(EDGE*id+RT+1) = 0.5*(pv_S+pv_LORT)
          dom%qe%elts(EDGE*id+DG+1) = 0.5*(pv_UPLT+pv_LORT)
          dom%qe%elts(EDGE*id+UP+1) = 0.5*(pv_UPLT+pv_W)
      end subroutine
  end subroutine

  function get_weights(dom, id, offs)
      type(Domain) dom
      real(8) get_weights(5)
      real(8) wgt(5)
      integer id, offs, i
      wgt(1) = dom%areas%elts(id+1)%part(1+offs)
      do i = 2, 5
          wgt(i) = wgt(i-1) + dom%areas%elts(id+1)%part(modulo(i+offs-1,6)+1)
      end do
      wgt = 0.5_8 - wgt*dom%areas%elts(id+1)%hex_inv
      get_weights = (/wgt(1), -wgt(2), wgt(3), -wgt(4), wgt(5)/)
  end function

  subroutine du_Qperp_Enstrophy(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer idNW
      integer idN
      integer idNE
      integer idW
      integer id
      integer idE
      integer idSW
      integer idS
      integer idSE
      real(8) wgt1(5), wgt2(5)
      idNW = idx(i - 1, j + 1, offs, dims)
      idN = idx(i, j + 1, offs, dims)
      idNE = idx(i + 1, j + 1, offs, dims)
      idW = idx(i - 1, j, offs, dims)
      id = idx(i, j, offs, dims)
      idE = idx(i + 1, j, offs, dims)
      idSW = idx(i - 1, j - 1, offs, dims)
      idS = idx(i, j - 1, offs, dims)
      idSE = idx(i + 1, j - 1, offs, dims)
      wgt1 = get_weights(dom, id, 0)
      wgt2 = get_weights(dom, idE, 3)
      dvelo(EDGE*id+RT+1) = dom%qe%elts(EDGE*id+RT+1)*( &
              sum(thickflux%data(dom%id+1)%elts((/ &
              DG+EDGE*id, EDGE*id+UP, EDGE*idW+RT, DG+EDGE*idSW, EDGE*idS+UP/)+1) * wgt1) + &
              sum(thickflux%data(dom%id+1)%elts((/ &
              DG+EDGE*idS, EDGE*idSE+UP, EDGE*idE+RT, DG+EDGE*idE, EDGE*idE+UP/)+1) * wgt2))
      wgt1 = get_weights(dom, id, 1)
      wgt2 = get_weights(dom, idNE, 4)
      dvelo(DG+EDGE*id+1) = dom%qe%elts(EDGE*id+DG+1)*( &
              sum(thickflux%data(dom%id+1)%elts((/ &
              EDGE*id+UP, EDGE*idW+RT, DG+EDGE*idSW, EDGE*idS+UP, EDGE*id+RT/)+1) * wgt1) + &
              sum(thickflux%data(dom%id+1)%elts((/ &
              EDGE*idE+UP, EDGE*idNE+RT, DG+EDGE*idNE, EDGE*idNE+UP, EDGE*idN+RT/)+1) * wgt2))
      wgt1 = get_weights(dom, id, 2)
      wgt2 = get_weights(dom, idN, 5)
      dvelo(EDGE*id+UP+1) = dom%qe%elts(EDGE*id+UP+1)*( &
              sum(thickflux%data(dom%id+1)%elts((/ &
              EDGE*idW+RT, DG+EDGE*idSW, EDGE*idS+UP, EDGE*id+RT, DG+EDGE*id/)+1) * wgt1) + &
              sum(thickflux%data(dom%id+1)%elts((/ &
              EDGE*idN+RT, DG+EDGE*idN, EDGE*idN+UP, EDGE*idNW+RT, DG+EDGE*idW/)+1) * wgt2))
  end subroutine

  subroutine du_source(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id, idE, idN, idNE
      real(8), dimension(RT:NODE) :: chi1, phi, hk_phi, kin_energ, u_mag
      id = idx(i, j, offs, dims)
      call du_Qperp(dom, i, j, offs, dims)
      if (viscosity .ne. 0) call diff_mom(dom, i, j, offs, dims)
      idN = idx(i, j + 1, offs, dims)
      idE = idx(i + 1, j, offs, dims)
      idNE = idx(i + 1, j + 1, offs, dims)
      if (penalize) then
          chi1 = penal%data(dom%id+1)%elts((/idN,idNE,idE,id/)+1)
          phi = 1 + alpha_m1*chi1
          chi1 = -ieta*chi1
          chi1(RT:UP) = 0.5_8*(chi1(NODE)+chi1(RT:UP)) ! interpolate p->u
          dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = dvelo(EDGE*id+RT+1:EDGE*id+UP+1) &
                  + velo(EDGE*id+RT+1:EDGE*id+UP+1)*dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1)*chi1(RT:UP)
      else
          phi = 1.0_8
      end if
      if (wind_stress .or. bottom_friction) then
          hk_phi = dom%topo%elts((/idN,idNE,idE,id/)+1)*phi + height((/idN,idNE,idE,id/)+1)
          where (hk_phi .ne. 0) hk_phi = phi/hk_phi
          hk_phi(RT:UP) = 0.5_8*(hk_phi(NODE)+hk_phi(RT:UP)) ! interpolate p->u
          if (wind_stress) then
              dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = dvelo(EDGE*id+RT+1:EDGE*id+UP+1) &
                      + dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1) &
                       *dom%windstress%elts(EDGE*id+RT+1:EDGE*id+UP+1)*hk_phi(RT:UP)
          end if
          if (bottom_friction) then
              kin_energ = dom%kin_energy%elts((/idN,idNE,idE,id/)+1)
              u_mag(RT:UP) = sqrt(2.0_8*0.5_8*(kin_energ(NODE)+kin_energ(RT:UP))) ! interpolate p->u
              dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = dvelo(EDGE*id+RT+1:EDGE*id+UP+1) &
                      - dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1) &
                       *friction_coeff*velo(EDGE*id+RT+1:EDGE*id+UP+1)*u_mag(RT:UP)*hk_phi(RT:UP)
          end if
      end if
  end subroutine

  subroutine du_Qperp(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer idNW
      integer idN
      integer idNE
      integer idW
      integer id
      integer idE
      integer idSW
      integer idS
      integer idSE
      real(8) wgt1(5), wgt2(5)
      idNW = idx(i - 1, j + 1, offs, dims)
      idN = idx(i, j + 1, offs, dims)
      idNE = idx(i + 1, j + 1, offs, dims)
      idW = idx(i - 1, j, offs, dims)
      id = idx(i, j, offs, dims)
      idE = idx(i + 1, j, offs, dims)
      idSW = idx(i - 1, j - 1, offs, dims)
      idS = idx(i, j - 1, offs, dims)
      idSE = idx(i + 1, j - 1, offs, dims)
      wgt1 = get_weights(dom, id, 0)
      wgt2 = get_weights(dom, idE, 3)
      dvelo(EDGE*id+RT+1) = &
              tflux(DG+EDGE*id+1)*0.5_8*(dom%qe%elts(DG+EDGE*id+1) + &
              dom%qe%elts(EDGE*id+RT+1))*wgt1(1) + &
              tflux(EDGE*id+UP+1)*0.5_8*(dom%qe%elts(EDGE*id+UP+1) + &
              dom%qe%elts(EDGE*id+RT+1))*wgt1(2) + &
              tflux(EDGE*idW+RT+1)*0.5_8*(dom%qe%elts(EDGE*idW+RT+1) + &
              dom%qe%elts(EDGE*id+RT+1))*wgt1(3) + &
              tflux(DG+EDGE*idSW+1)*0.5_8*(dom%qe%elts(DG+EDGE*idSW+1) &
              + dom%qe%elts(EDGE*id+RT+1))*wgt1(4) + &
              tflux(EDGE*idS+UP+1)*0.5_8*(dom%qe%elts(EDGE*idS+UP+1) + &
              dom%qe%elts(EDGE*id+RT+1))*wgt1(5) + &
              tflux(DG+EDGE*idS+1)*0.5_8*(dom%qe%elts(DG+EDGE*idS+1) + &
              dom%qe%elts(EDGE*id+RT+1))*wgt2(1) + &
              tflux(EDGE*idSE+UP+1)*0.5_8*(dom%qe%elts(EDGE*idSE+UP+1) &
              + dom%qe%elts(EDGE*id+RT+1))*wgt2(2) + &
              tflux(EDGE*idE+RT+1)*0.5_8*(dom%qe%elts(EDGE*idE+RT+1) + &
              dom%qe%elts(EDGE*id+RT+1))*wgt2(3) + &
              tflux(DG+EDGE*idE+1)*0.5_8*(dom%qe%elts(DG+EDGE*idE+1) + &
              dom%qe%elts(EDGE*id+RT+1))*wgt2(4) + &
              tflux(EDGE*idE+UP+1)*0.5_8*(dom%qe%elts(EDGE*idE+UP+1) + &
              dom%qe%elts(EDGE*id+RT+1))*wgt2(5)
      wgt1 = get_weights(dom, id, 1)
      wgt2 = get_weights(dom, idNE, 4)
      dvelo(DG+EDGE*id+1) = &
              tflux(EDGE*id+UP+1)*0.5_8*(dom%qe%elts(EDGE*id+UP+1) + &
              dom%qe%elts(DG+EDGE*id+1))*wgt1(1) + &
              tflux(EDGE*idW+RT+1)*0.5_8*(dom%qe%elts(EDGE*idW+RT+1) + &
              dom%qe%elts(DG+EDGE*id+1))*wgt1(2) + &
              tflux(DG+EDGE*idSW+1)*0.5_8*(dom%qe%elts(DG+EDGE*idSW+1) &
              + dom%qe%elts(DG+EDGE*id+1))*wgt1(3) + &
              tflux(EDGE*idS+UP+1)*0.5_8*(dom%qe%elts(EDGE*idS+UP+1) + &
              dom%qe%elts(DG+EDGE*id+1))*wgt1(4) + &
              tflux(EDGE*id+RT+1)*0.5_8*(dom%qe%elts(EDGE*id+RT+1) + &
              dom%qe%elts(DG+EDGE*id+1))*wgt1(5) + &
              tflux(EDGE*idE+UP+1)*0.5_8*(dom%qe%elts(EDGE*idE+UP+1) + &
              dom%qe%elts(DG+EDGE*id+1))*wgt2(1) + &
              tflux(EDGE*idNE+RT+1)*0.5_8*(dom%qe%elts(EDGE*idNE+RT+1) &
              + dom%qe%elts(DG+EDGE*id+1))*wgt2(2) + &
              tflux(DG+EDGE*idNE+1)*0.5_8*(dom%qe%elts(DG+EDGE*idNE+1) &
              + dom%qe%elts(DG+EDGE*id+1))*wgt2(3) + &
              tflux(EDGE*idNE+UP+1)*0.5_8*(dom%qe%elts(EDGE*idNE+UP+1) &
              + dom%qe%elts(DG+EDGE*id+1))*wgt2(4) + &
              tflux(EDGE*idN+RT+1)*0.5_8*(dom%qe%elts(EDGE*idN+RT+1) + &
              dom%qe%elts(DG+EDGE*id+1))*wgt2(5)
      wgt1 = get_weights(dom, id, 2)
      wgt2 = get_weights(dom, idN, 5)
      dvelo(EDGE*id+UP+1) = &
              tflux(EDGE*idW+RT+1)*0.5_8*(dom%qe%elts(EDGE*idW+RT+1) + &
              dom%qe%elts(EDGE*id+UP+1))*wgt1(1) + &
              tflux(DG+EDGE*idSW+1)*0.5_8*(dom%qe%elts(DG+EDGE*idSW+1) &
              + dom%qe%elts(EDGE*id+UP+1))*wgt1(2) + &
              tflux(EDGE*idS+UP+1)*0.5_8*(dom%qe%elts(EDGE*idS+UP+1) + &
              dom%qe%elts(EDGE*id+UP+1))*wgt1(3) + &
              tflux(EDGE*id+RT+1)*0.5_8*(dom%qe%elts(EDGE*id+RT+1) + &
              dom%qe%elts(EDGE*id+UP+1))*wgt1(4) + &
              tflux(DG+EDGE*id+1)*0.5_8*(dom%qe%elts(DG+EDGE*id+1) + &
              dom%qe%elts(EDGE*id+UP+1))*wgt1(5) + &
              tflux(EDGE*idN+RT+1)*0.5_8*(dom%qe%elts(EDGE*idN+RT+1) + &
              dom%qe%elts(EDGE*id+UP+1))*wgt2(1) + &
              tflux(DG+EDGE*idN+1)*0.5_8*(dom%qe%elts(DG+EDGE*idN+1) + &
              dom%qe%elts(EDGE*id+UP+1))*wgt2(2) + &
              tflux(EDGE*idN+UP+1)*0.5_8*(dom%qe%elts(EDGE*idN+UP+1) + &
              dom%qe%elts(EDGE*id+UP+1))*wgt2(3) + &
              tflux(EDGE*idNW+RT+1)*0.5_8*(dom%qe%elts(EDGE*idNW+RT+1) &
              + dom%qe%elts(EDGE*id+UP+1))*wgt2(4) + &
              tflux(DG+EDGE*idW+1)*0.5_8*(dom%qe%elts(DG+EDGE*idW+1) + &
              dom%qe%elts(EDGE*id+UP+1))*wgt2(5)
  end subroutine

  subroutine height_trend(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id
      integer idS
      integer idW
      integer idSW
      id = idx(i, j, offs, dims)
      idS = idx(i, j - 1, offs, dims)
      idW = idx(i - 1, j, offs, dims)
      idSW = idx(i - 1, j - 1, offs, dims)
      dheight(id+1) = -(tflux(EDGE*id+UP+1)  - tflux(DG+EDGE*id+1)   + tflux(EDGE*id+RT+1) - &
                       tflux(EDGE*idS+UP+1) + tflux(DG+EDGE*idSW+1) - tflux(EDGE*idW+RT+1)) &
              *dom%areas%elts(id+1)%hex_inv
  end subroutine

  subroutine du_gradB(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id
      integer idE
      integer idN
      integer idNE
      id = idx(i, j, offs, dims)
      idE = idx(i + 1, j, offs, dims)
      idN = idx(i, j + 1, offs, dims)
      idNE = idx(i + 1, j + 1, offs, dims)
      dvelo(EDGE*id+RT+1) = (dvelo(EDGE*id+RT+1) - &
              (dom%bernoulli%elts(idE+1) - &
              dom%bernoulli%elts(id+1)))/dom%len%elts(EDGE*id+RT+1)
      dvelo(DG+EDGE*id+1) = (dvelo(DG+EDGE*id+1) - &
              (dom%bernoulli%elts(id+1) - &
              dom%bernoulli%elts(idNE+1)))/dom%len%elts(DG+EDGE*id+1)
      dvelo(EDGE*id+UP+1) = (dvelo(EDGE*id+UP+1) - &
              (dom%bernoulli%elts(idN+1) - &
              dom%bernoulli%elts(id+1)))/dom%len%elts(EDGE*id+UP+1)
  end subroutine

  subroutine comp_offs3(dom, p, offs, dims)
    type(Domain) dom
    integer p 
    integer offs(0:N_BDRY), dims(2,N_BDRY), i, n
    offs(0) = dom%patch%elts(p+1)%elts_start
    do i = 1, N_BDRY
        n = dom%patch%elts(p+1)%neigh(i)
        if (n .gt. 0) then ! regular patch
            offs(i)  = dom%patch%elts(n+1)%elts_start
            dims(:,i) = PATCH_SIZE
        else if (n .lt. 0) then
            offs(i)  = dom%bdry_patch%elts(-n+1)%elts_start
            dims(:,i) = sides_dims(:,abs(dom%bdry_patch%elts(-n+1)%side)+1)
        end if
    end do
    offs(1:N_BDRY) = offs(1:N_BDRY) - offs(0) ! make relative
    offs(SOUTH) = offs(SOUTH) +  dims(1,SOUTH)*(dims(2,SOUTH)-1)
    offs(SOUTHEAST) = offs(SOUTHEAST) + dims(1,SOUTHEAST)*(dims(2,SOUTH)-1)
    offs(WEST)  = offs(WEST) + (dims(1,WEST)-1)
    offs(NORTHWEST) = offs(NORTHWEST) + (dims(1,NORTHWEST)-1)
    offs(SOUTHWEST) = offs(SOUTHWEST) + dims(1,SOUTHWEST)*dims(2,SOUTHWEST) - 1
    offs(NORTH) = offs(NORTH) - PATCH_SIZE*(PATCH_SIZE-1) 
    offs(NORTHWEST) = offs(NORTHWEST) - PATCH_SIZE*(PATCH_SIZE-1) 
    offs(EAST)  = offs(EAST)  -            (PATCH_SIZE-1)
    offs(SOUTHEAST) = offs(SOUTHEAST) -(PATCH_SIZE-1)
    offs(NORTHEAST) = offs(NORTHEAST) - (PATCH_SIZE*PATCH_SIZE-1)
  end subroutine
end module ops_mod
