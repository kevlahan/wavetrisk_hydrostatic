"""Stage226 production integration and bounded higher-ratio correctness gates."""
from capability_223 import digest, save, verify_kit, verify_sources, coverage, events, numerical


def schedule():
    # Anchor the production integration at the already validated boundary.
    yield ('J4',40,'baseline',False,'native')
    yield ('J4',40,'candidate',False,'native')
    for n in (41,80,120,160):
        yield ('J4',n,'checked',False,'native')
        yield ('J4',n,'candidate',False,'native')
        if n in (41,160):yield ('J4',n,'candidate',True,'native')
        if n==41:yield ('J4',n,'checked',False,'empty-root')
    yield ('J5',160,'baseline',False,'native')
    yield ('J5',160,'candidate',False,'native')
    for n in (167,256):
        yield ('J5',n,'checked',False,'native')
        yield ('J5',n,'candidate',False,'native')
        yield ('J5',n,'candidate',True,'native')
