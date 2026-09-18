#!/usr/bin/env python3
"""Remove unused scalar work from the two velocity-only reconstruction sites.

Applied only to the frozen, uninstrumented Stage202 tree. No arithmetic change,
no Domain-authority shortcut, and no persistent plan/epoch change.
"""
from pathlib import Path


def once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Source contract changed: '+old[:100])
    return text.replace(old, new)


def candidate(root):
    path = root/'src/parallel_block_mpi.f90'
    source = path.read_text()
    for start in ('validate_candidate_block_velocity_restriction',
                  'validate_complete_block_vector_wavelets'):
        lo = source.index('  subroutine '+start+' (')
        hi = source.index('  end subroutine '+start, lo)
        part = source[lo:hi]
        part = once(part, 'call reconstruct_block_writeback_domain_stage(BLOCK_PAYLOAD_SOL)',
                    'call reconstruct_block_writeback_domain_stage( &\n            BLOCK_PAYLOAD_SOL,vector_only=.true.)')
        if 'scalar_domain_stage' in part:
            raise ValueError('Vector consumer unexpectedly reads scalar staging')
        source = source[:lo]+part+source[hi:]
    lo = source.index('  subroutine exchange_block_writeback_payloads (')
    hi = source.index('  end subroutine exchange_block_writeback_payloads', lo)
    part = source[lo:hi]
    part = once(part, 'payload_family,component_family)', 'payload_family,component_family,skip_scalar_clear)')
    part = once(part, 'integer, optional, intent(in) :: component_family',
                'integer, optional, intent(in) :: component_family\n    logical, optional, intent(in) :: skip_scalar_clear\n    logical :: omit_scalar_clear')
    part = once(part, '    block_writeback_plan%scalar_send_buffer = 0.0_dp\n    block_writeback_plan%scalar_recv_buffer = 0.0_dp',
                '''    omit_scalar_clear = .false.
    if (present(skip_scalar_clear)) omit_scalar_clear = skip_scalar_clear
    if (omit_scalar_clear .and. component /= BLOCK_WRITEBACK_VECTOR) &
         call fail("Stage209 scalar omission requires vector-only exchange")
    if (.not. omit_scalar_clear) then
       block_writeback_plan%scalar_send_buffer = 0.0_dp
       block_writeback_plan%scalar_recv_buffer = 0.0_dp
    end if''')
    source = source[:lo]+part+source[hi:]
    lo = source.index('  subroutine reconstruct_block_writeback_domain_stage (')
    hi = source.index('  end subroutine reconstruct_block_writeback_domain_stage', lo)
    part = source[lo:hi]
    part = once(part, '  subroutine reconstruct_block_writeback_domain_stage ( &\n       payload_family,component_family)', '  subroutine reconstruct_block_writeback_domain_stage ( &\n       payload_family,component_family,vector_only)')
    part = once(part, 'integer, optional, intent(in) :: component_family',
                '''integer, optional, intent(in) :: component_family
    logical, optional, intent(in) :: vector_only
    logical :: omit_scalar
    integer :: component''')
    part = once(part, '''    call exchange_block_writeback_payloads( &
         payload_family,component_family)

    block_writeback_plan%scalar_domain_stage = 0.0_dp''', '''    omit_scalar = .false.
    if (present(vector_only)) omit_scalar = vector_only
    component = BLOCK_WRITEBACK_BOTH
    if (present(component_family)) component = component_family
    if (omit_scalar) then
       if (present(component_family)) &
            call fail("Stage209 vector-only reconstruction has two selectors")
       component = BLOCK_WRITEBACK_VECTOR
    end if
    call exchange_block_writeback_payloads( &
         payload_family,component,skip_scalar_clear=omit_scalar)

    if (.not. omit_scalar) block_writeback_plan%scalar_domain_stage = 0.0_dp
#ifdef STAGE209_POISON
    ! Validation-only cut: omitted arrays must never supply a consumer.
    if (omit_scalar) then
       block_writeback_plan%scalar_domain_stage = -huge(0.0_dp)
       block_writeback_plan%scalar_send_buffer = -huge(0.0_dp)
       block_writeback_plan%scalar_recv_buffer = -huge(0.0_dp)
    end if
#endif''')
    part = once(part, '''                     vector_stage_start)
                block_writeback_plan%domain_patch_covered''',
                '''                     vector_stage_start,skip_scalar=omit_scalar)
                block_writeback_plan%domain_patch_covered''')
    part = once(part, '      call get_local_block_scalar_patch_family_values( &',
                '      if (.not. omit_scalar) call get_local_block_scalar_patch_family_values( &')
    part = once(part, '''      block_writeback_plan%scalar_domain_stage( &
           scalar_stage_start:scalar_stage_start+ &
           block_writeback_plan%scalar_patch_nvalue-1) = &
           scalar_payload(scalar_pos:scalar_pos+ &
           block_writeback_plan%scalar_patch_nvalue-1)''',
                '''      if (.not. omit_scalar) then
         block_writeback_plan%scalar_domain_stage( &
              scalar_stage_start:scalar_stage_start+ &
              block_writeback_plan%scalar_patch_nvalue-1) = &
              scalar_payload(scalar_pos:scalar_pos+ &
              block_writeback_plan%scalar_patch_nvalue-1)
      end if''')
    source = source[:lo]+part+source[hi:]
    source = once(source,
        '    if (omit_scalar.and.payload_family/=BLOCK_PAYLOAD_WAV_COEFF) &\n         call fail("scalar-free packing requires wavelet coefficients")',
        '    if (omit_scalar .and. payload_family /= BLOCK_PAYLOAD_WAV_COEFF .and. &\n         payload_family /= BLOCK_PAYLOAD_SOL) &\n         call fail("scalar-free packing requires SOL or wavelet coefficients")')
    path.write_text(source)

if __name__ == '__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('tree',type=Path)
    candidate(p.parse_args().tree)
