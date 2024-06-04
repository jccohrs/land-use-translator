import unittest
from src.lut import LUT

class TestClass(unittest.TestCase):
    def test_namelist(self):
        config = dotdict(main())

        LUT = LUT(config)
        namelist = LUT.generate_namelist()
        print(namelist)
        assert namelist is True

if __name__=='__main__':
	unittest.main()
