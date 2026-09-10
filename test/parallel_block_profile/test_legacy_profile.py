import unittest
from prepare_legacy_profile import replace_once, scope, wrap_calls


class TransformTests(unittest.TestCase):
    def test_no_argument_call_does_not_consume_next_line(self):
        text = '    call restart\n\n  end subroutine write_checkpoint\n'
        wrapped = wrap_calls(text, 'restart', 'DP_RESTART')
        self.assertIn('call restart\n    call detail_leave(DP_RESTART)\n\n  end', wrapped)

    def test_call_name_is_exact(self):
        text = '    call RK_sub_step (q, h)\n    call RK_sub_step2 (q, h)\n'
        wrapped = wrap_calls(text, 'RK_sub_step', 'DP_RK_ASSEMBLE')
        self.assertEqual(wrapped.count('call detail_enter'), 1)

    def test_scope_rejects_return(self):
        with self.assertRaises(ValueError):
            scope('  subroutine f\n    if (x) return\n  end subroutine f', 'f', '    if', 'DP_STEP')

    def test_anchor_requires_unique_match(self):
        with self.assertRaises(ValueError):
            replace_once('xx', 'x', 'y')


if __name__ == '__main__':
    unittest.main()
