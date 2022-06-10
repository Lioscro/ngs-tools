import shutil
from unittest import TestCase

from ngs_tools import binary

from . import mixins


class TestBinary(mixins.TestMixin, TestCase):

    def test_ls(self):
        bin = binary.Binary('ls')
        executor = binary.BinaryExecutor(
            bin,
            binary.PositionalArgument('file'),
            binary.ConstantArgument('-a', required=False),
        )
        result = executor({'-a': None, 'file': '.'})
        self.assertEqual(result.command, [shutil.which('ls'), '-a', '.'])

        result = executor({'file': '.'})
        self.assertEqual(result.command, [shutil.which('ls'), '.'])

        with self.assertRaises(binary.BinaryError):
            executor({'file': '.', 'other': None})

        with self.assertRaises(binary.BinaryError):
            executor = binary.BinaryExecutor(
                bin,
                binary.PositionalArgument('file'),
                binary.ConstantArgument('-a', required=True),
            )
            executor({'file': '.'})


class TestArgumentValidator(mixins.TestMixin, TestCase):

    def test_str(self):
        self.assertEqual(
            str(binary.IsPositiveInteger), '(IsInteger AND IsPositive)'
        )

    def test_int(self):
        self.assertTrue(binary.IsInteger('3'))
        self.assertFalse(binary.IsInteger('3.2'))

    def test_float(self):
        self.assertTrue(binary.IsFloat('3.2'))
        self.assertFalse(binary.IsFloat('a'))

    def test_positive(self):
        self.assertTrue(binary.IsPositive('3'))
        self.assertFalse(binary.IsPositive('-3.2'))

    def test_positive_integer(self):
        self.assertTrue(binary.IsPositiveInteger('3'))
        self.assertFalse(binary.IsPositiveInteger('-3'))
        self.assertFalse(binary.IsPositiveInteger('3.2'))
