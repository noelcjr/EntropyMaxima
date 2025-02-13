import unittest, codecs
from em.charmm.gen import gen_setup_one

class TestGen(unittest.TestCase):
    def test_numbers_3_4(self):
        with codecs.open('../../../charmm_templates/setup_one.inp','r', encoding='utf-8') as f:
            expected = f.read()
            result = gen_setup_one()
            self.maxDiff = None
            self.assertEqual(expected, result)

if __name__ == '__main__':
    unittest.main()
