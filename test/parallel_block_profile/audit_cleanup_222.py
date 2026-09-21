"""Add bounded numerical traces and deletion poison to disposable source trees.

No timing instrumentation or solver-side audit switches are added to production.
Reference and candidate emit the first RK4 step's Domain/compact transitions.
The candidate additionally poisons omitted inverse stores and unused staging on
every step. Run from a fresh directory: stream files are append-only.
"""
import argparse
from pathlib import Path
import re


def once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f'Expected one anchor: {old[:100]} ({text.count(old)})')
    return text.replace(old, new)


def routine(text, name, transform):
    start = re.search(r'^  subroutine '+name+r'\b', text, re.M)
    end = re.search(r'^  end subroutine '+name+r'\b', text, re.M)
    if not start or not end:
        raise ValueError(name)
    return text[:start.start()]+transform(text[start.start():end.end()])+text[end.end():]


TRACE = '''
  subroutine audit_values222(label,values)
    character(*),intent(in) :: label
    real(dp),intent(in) :: values(:)
    integer :: unit
    character(128) :: filename
    write(filename,'(a,a,a,i0,a)') 'cleanup222-',trim(label),'-',rank,'.bin'
    open(newunit=unit,file=trim(filename),access='stream',form='unformatted', &
         status='unknown',position='append')
    write(unit) size(values),values
    close(unit)
  end subroutine audit_values222

  subroutine audit_field222(label,fields)
    character(*),intent(in) :: label
    type(Float_Field),intent(in) :: fields(1:N_VARIABLE,1:zlevels)
    integer,save :: calls(3)=0
    integer :: which,d,k,v,n
    select case(label)
    case('inverse');which=1
    case('rk');which=2
    case('forward');which=3
    case default;error stop 'unexpected cleanup audit field'
    end select
    calls(which)=calls(which)+1
    if(calls(which)>4) return
    do d=1,size(grid)
       do k=1,zlevels
          do v=1,N_VARIABLE
             n=fields(v,k)%data(d)%length
             call audit_values222(label,fields(v,k)%data(d)%elts(1:n))
          end do
       end do
    end do
  end subroutine audit_field222

  subroutine audit_blocks222(family)
    integer,intent(in) :: family
    integer,save :: calls=0
    integer :: slot,b,p,n
    real(dp),allocatable :: values(:)
    calls=calls+1
    if(calls>4) return
    do slot=1,n_local_blocks()
       b=local_block_catalog(slot)
       n=local_block_scalar_family_patch_nvalue(b)
       allocate(values(n))
       do p=0,local_block_patch_count(b)-1
          call get_local_block_scalar_patch_family_values(b,p,family,values)
          call audit_values222('compact_scalar',values)
       end do
       deallocate(values)
       n=local_block_vector_family_patch_nvalue(b)
       allocate(values(n))
       do p=0,local_block_patch_count(b)-1
          call get_local_block_vector_patch_family_values(b,p,family,values)
          call audit_values222('compact_vector',values)
       end do
       deallocate(values)
    end do
  end subroutine audit_blocks222
'''


def instrument(root, poison=False):
    def change(name, transform):
        path = root/'src'/name
        path.write_text(transform(path.read_text()))

    def mpi(s):
        s = s.replace('\n  implicit none\n',
                      '\n  implicit none\n  public :: audit_values222, audit_field222\n', 1)
        end = s.lower().rfind('end module')
        s = s[:end]+TRACE+s[end:]
        s = routine(s, 'activate_block_native_inverse_transform', lambda t: once(t,
                    '    if (validate_oracle) then\n       call block_profile_enter(BLOCK_PROFILE_ORACLE)\n       call assert_block_domain_field_family_match(',
                    "    call audit_field222('inverse',domain_scaling)\n"+
                    '    if (validate_oracle) then\n       call block_profile_enter(BLOCK_PROFILE_ORACLE)\n       call assert_block_domain_field_family_match('))
        s = routine(s, 'retain_block_native_multistage_candidate', lambda t: once(t,
                    '    call detail_leave(DP_RETAIN)',
                    "    call audit_field222('rk',domain_sol)\n    call audit_blocks222(BLOCK_PAYLOAD_SOL)\n    call detail_leave(DP_RETAIN)"))
        s = routine(s, 'activate_block_native_wavelet_compression', lambda t: once(t,
                    '  end subroutine activate_block_native_wavelet_compression',
                    "    call audit_field222('forward',domain_wavelet)\n  end subroutine activate_block_native_wavelet_compression"))
        if poison:
            s = once(s, '  use ieee_arithmetic, only : ieee_is_finite',
                     '  use ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan')
            s = once(s, '  use parallel_block_inverse_mod, only : prepare_native_inverse',
                     '  use parallel_block_inverse_mod, only : poison_publication222\n  use parallel_block_inverse_mod, only : prepare_native_inverse')
            s = once(s, '    if (.not.provisional_transform.or.validate_oracle) &',
                     '    if (provisional_transform.and..not.validate_oracle) call poison_publication222(domain_scaling)\n'+
                     '    if (.not.provisional_transform.or.validate_oracle) &')
            s = routine(s, 'reconstruct_block_writeback_domain_stage', lambda t: once(t,
                        '\n    if (temperature_only) then',
                        '''
    if (temperature_only) then
       block_writeback_plan%scalar_domain_stage = ieee_value(0.0_dp,ieee_quiet_nan)
       block_writeback_plan%vector_domain_stage = ieee_value(0.0_dp,ieee_quiet_nan)
       write(*,'(a,i0)') 'Stage222 poisoned staging rank = ',rank'''))
            s = routine(s, 'prepare_parallel_block_grid_change', lambda t: once(t,
                        '    if (.not.native_publication_ready) call synchronize_parallel_block_checkpoint(.true.)',
                        '''    if(native_publication_ready) then
       call assert_block_domain_field_family_match(BLOCK_PAYLOAD_SOL)
       call assert_block_domain_field_family_match(BLOCK_PAYLOAD_WAV_COEFF)
       write(*,'(a,i0)') 'Stage222 skipped grid publication verified rank = ',rank
    end if
    if (.not.native_publication_ready) call synchronize_parallel_block_checkpoint(.true.)'''))
        return s
    change('parallel_block_mpi.f90', mpi)

    def time(s):
        s = s.replace('\n', '\n  use parallel_block_mpi_mod, only : audit_values222\n', 1)
        def tendency(t):
            t = once(t, '    real(dp) :: profile_start',
                     '    real(dp) :: profile_start\n    integer,save :: audit_calls222=0\n    integer :: d222,k222')
            return once(t, '    call begin_block_domain_multistage_candidate_stage( &',
                        '''    audit_calls222=audit_calls222+1
    if(audit_calls222<=4) then
       do d222=1,size(grid)
          do k222=1,zlevels
             call audit_values222('mass_tendency',native_mass(d222)%tendency(:,k222))
             call audit_values222('velocity_tendency',native_velocity(d222)%tendency(:,k222))
          end do
       end do
    end if
    call begin_block_domain_multistage_candidate_stage( &''')
        return routine(s, 'prepare_block_multistage_tendency', tendency)
    change('time_integr.f90', time)

    if poison:
        def inverse(s):
            s = once(s, '  use ieee_arithmetic, only : ieee_is_finite',
                     '  use ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan')
            s = once(s, '  public :: capture_native_inverse_consumer_aliases',
                     '  public :: poison_publication222\n  public :: capture_native_inverse_consumer_aliases')
            s = once(s, '\ncontains\n', '''
contains
  subroutine poison_publication222(scaling)
    type(Float_Field),intent(inout) :: scaling(1:N_VARIABLE,1:zlevels)
    integer :: p,d,id,k,v,n
    n=0
    do p=1,nnode
       d=node_domain(p);id=node_id(p)
       if(nodes(d)%key(1,id+1)>0) cycle
       n=n+1
       do k=1,zlevels
          do v=1,nscalar
             scaling(scalars(1)+v-1,k)%data(d)%elts(id+1)=ieee_value(0.0_dp,ieee_quiet_nan)
          end do
          do v=1,EDGE
             scaling(S_VELO,k)%data(d)%elts(EDGE*id+v)=ieee_value(0.0_dp,ieee_quiet_nan)
          end do
       end do
    end do
    write(*,'(a,i0,a,i0)') 'Stage222 poisoned inverse rank = ',rank,' slots = ',n
  end subroutine poison_publication222
''')
            return s
        change('parallel_block_inverse.f90', inverse)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tree', type=Path)
    parser.add_argument('--poison', action='store_true')
    args = parser.parse_args()
    instrument(args.tree, args.poison)
