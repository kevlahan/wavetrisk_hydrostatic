import unittest
from migrate_scalar_storage_access import migrate,transform


class StorageMigrationTests(unittest.TestCase):
    def test_scalar_range_vector_and_inline_if(self):
        ref='block_scalar_tendency(i)%patch'
        self.assertEqual(transform(ref+'(p)=a'),f'call scalar_write({ref},p,a)')
        self.assertEqual(transform(ref+'(p:q)=a'),f'call scalar_write_range({ref},p,q,a)')
        self.assertEqual(transform('if (ready) '+ref+'(p:q)=a'),f'if (ready) call scalar_write_range({ref},p,q,a)')
        self.assertEqual(transform('v='+ref+'(base+index)'),f'v=scalar_read({ref},base+index)')

    def test_preserve_other_code_and_idempotence(self):
        text='  ! unchanged comment\n  call scalar_fill(block_scalar_tendency(i)%patch,0.0_dp)\n'
        self.assertEqual(migrate(text),(text,0))
        text='  x=block_scalar_tendency(i)%patch( &\n       a:b)\n'
        result,n=migrate(text)
        self.assertEqual(n,1)
        self.assertEqual(migrate(result),(result,0))


if __name__=='__main__':unittest.main()
