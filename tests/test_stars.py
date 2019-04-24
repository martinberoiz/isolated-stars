import unittest
import stars

class TestAlign(unittest.TestCase):
    def test_isolated(self):
        stars.isolated([], 10.0)

if __name__ == "__main__":
    unittest.main()
