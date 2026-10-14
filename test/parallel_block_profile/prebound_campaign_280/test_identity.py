import json,unittest
from pathlib import Path
from common import require_identity,REVISION
import run,review
class Identity(unittest.TestCase):
    def test_reader_and_runner_share_sealed_revision(self):
        identity=json.loads((Path(__file__).parent/'IDENTITY.json').read_text())
        require_identity(identity)
        self.assertEqual(run.REVISION,REVISION);self.assertEqual(review.REVISION,REVISION)
        for key,value in [('revision','resident-transform-v1'),('stage',278),('main_commit','changed')]:
            with self.subTest(key=key),self.assertRaises(ValueError):require_identity(dict(identity,**{key:value}))
