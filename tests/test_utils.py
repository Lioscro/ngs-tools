import gzip
import os
import pickle
import subprocess
from unittest import TestCase, mock

from joblib import delayed

from ngs_tools import utils

from . import mixins


class TestUtils(mixins.TestMixin, TestCase):

    def test_retry(self):
        test = mock.MagicMock()
        test.return_value = 'test'
        self.assertEqual(
            'test',
            utils.retry(test, retries=3, args=('a', 'b'), kwargs={'c': 'd'})
        )
        test.assert_called_once_with('a', 'b', c='d')

    def test_retry_raises(self):
        test = mock.MagicMock()
        test.side_effect = Exception('test')
        with self.assertRaises(Exception):
            utils.retry(test, retries=3)
        self.assertEqual(3, test.call_count)

    def test_retry_decorator(self):

        def test_func(*args, **kwargs):
            return

        with mock.patch('ngs_tools.utils.retry') as retry:
            decorated_func = utils.retry_decorator(3, 1, True, [Exception])(
                test_func
            )
            decorated_func('t', kw='arg')

            retry.assert_called_once_with(
                test_func,
                3, ('t',), {'kw': 'arg'},
                retry_every=1,
                backoff=True,
                exceptions=[Exception]
            )

    @utils.retry_decorator(3)
    def test_run_executable(self):
        p, stdout, stderr = utils.run_executable(['echo', 'TEST'], stream=False)
        self.assertEqual('TEST\n', stdout)

    def test_run_exectuable_raises_exception(self):
        with self.assertRaises(subprocess.SubprocessError):
            utils.run_executable(['bash', 'nonexistent option'])

    def test_run_exectuable_with_returncode(self):
        utils.run_executable(['bash', 'nonexistent option'], returncode=127)

    def test_run_executable_no_wait(self):
        with mock.patch('ngs_tools.utils.subprocess') as sp_mock:
            sp_mock.Popen().returncode = 0
            utils.run_executable(['echo', 'TEST'], wait=False)
            sp_mock.Popen().poll.assert_not_called()

    @utils.retry_decorator(3)
    def test_run_executable_with_stream(self):
        with mock.patch('ngs_tools.utils.logger.debug') as debug_mock:
            utils.run_executable(['echo', 'TEST'], stream=True)
            debug_mock.assert_has_calls([mock.call('TEST')])

    def test_ParallelWithProgress(self):
        utils.ParallelWithProgress(
            delayed(mixins.dummy_function)() for _ in range(10)
        )

    def test_is_gzip(self):
        self.assertTrue(utils.is_gzip(self.fastq_gz_path))
        self.assertFalse(utils.is_gzip(self.fastq_path))

    def test_is_remote(self):
        self.assertFalse(utils.is_remote('path/to/local/file'))
        self.assertTrue(utils.is_remote('https://path/to/remote/file'))

    def test_mkstemp(self):
        path = utils.mkstemp()
        self.assertTrue(os.path.exists(path))

    def test_mkstemp_delete(self):
        path = utils.mkstemp(delete=True)
        self.assertFalse(os.path.exists(path))

    def test_decompress_gzip(self):
        gzip_path = os.path.join(self.temp_dir, 'compressed.gz')
        out_path = os.path.join(self.temp_dir, 'decompressed')
        with gzip.open(gzip_path, 'wt') as f:
            f.write('TESTING\nTEST')
        self.assertEqual(out_path, utils.decompress_gzip(gzip_path, out_path))
        self.assertTrue(os.path.exists(out_path))
        with open(out_path, 'r') as f:
            self.assertEqual('TESTING\nTEST', f.read())

    def test_compress_gzip(self):
        file_path = os.path.join(self.temp_dir, 'file')
        out_path = os.path.join(f'{file_path}.gz')
        with open(file_path, 'w') as f:
            f.write('TESTING\nTEST')
        self.assertEqual(out_path, utils.compress_gzip(file_path, out_path))
        self.assertTrue(os.path.exists(out_path))
        with gzip.open(out_path, 'rt') as f:
            self.assertEqual('TESTING\nTEST', f.read())

    def test_concatenate_files(self):
        file1_path = os.path.join(self.temp_dir, 'file1')
        file2_path = os.path.join(self.temp_dir, 'file2')

        with open(file1_path, 'w') as f:
            f.write('TEST1\n')
        with open(file2_path, 'w') as f:
            f.write('TEST2\n')

        out_path = utils.concatenate_files(
            file1_path,
            file2_path,
            out_path=os.path.join(self.temp_dir, 'concatenated'),
        )

        with open(out_path, 'r') as f:
            self.assertEqual(f.read(), 'TEST1\nTEST2\n')

    def test_concatenate_files_as_text(self):
        file1_path = os.path.join(self.temp_dir, 'file1')
        file2_path = os.path.join(self.temp_dir, 'file2.gz')

        with open(file1_path, 'w') as f:
            f.write('TEST1')
        with gzip.open(file2_path, 'wt') as f:
            f.write('TEST2')

        out_path = utils.concatenate_files_as_text(
            file1_path,
            file2_path,
            out_path=os.path.join(self.temp_dir, 'concatenated'),
        )

        with open(out_path, 'r') as f:
            self.assertEqual(f.read(), 'TEST1\nTEST2\n')

    def test_download_file(self):
        with mock.patch('ngs_tools.utils.urlretrieve') as urlretrieve,\
            mock.patch('ngs_tools.utils.TqdmUpTo', mixins.tqdm_mock):
            path = os.path.join(self.temp_dir, 'test')
            self.assertEqual(path, utils.download_file('remote/path', path))
            urlretrieve.assert_called_once_with(
                'remote/path',
                filename=path,
                reporthook=mock.ANY,
                data=mock.ANY
            )

    def test_stream_file(self):
        with mock.patch('ngs_tools.utils.os.mkfifo'),\
            mock.patch('ngs_tools.utils.urlretrieve') as urlretrieve,\
            mock.patch('ngs_tools.utils.threading.Thread') as Thread:
            with utils.stream_file('remote/path', 'local/path') as path:
                self.assertEqual('local/path', path)
                Thread.assert_called_once_with(
                    target=urlretrieve,
                    args=('remote/path', 'local/path'),
                    daemon=True
                )
                Thread.return_value.start.assert_called_once()
                Thread.return_value.join.assert_not_called()
            Thread.return_value.join.assert_called_once()

    def test_stream_file_not_supported(self):
        with mock.patch('ngs_tools.utils.os.mkfifo') as mkfifo,\
            mock.patch('ngs_tools.utils.urlretrieve'),\
            mock.patch('ngs_tools.utils.threading.Thread') as Thread:
            mkfifo.side_effect = AttributeError('test')
            with self.assertRaises(OSError):
                with utils.stream_file('remote/path', 'local/path'):
                    pass
                Thread.assert_not_called()

    def test_stream_file_exception(self):
        with mock.patch('ngs_tools.utils.os.mkfifo'),\
            mock.patch('ngs_tools.utils.urlretrieve') as urlretrieve,\
            mock.patch('ngs_tools.utils.threading.Thread') as Thread:
            with self.assertRaises(Exception):
                with utils.stream_file('remote/path', 'local/path') as path:
                    self.assertEqual('local/path', path)
                    Thread.assert_called_once_with(
                        target=urlretrieve,
                        args=('remote/path', 'local/path'),
                        daemon=True
                    )
                    Thread.return_value.start.assert_called_once()
                    Thread.return_value.join.assert_not_called()
                    raise Exception('test')
            Thread.return_value.join.assert_called_once()

    def test_write_pickle(self):
        path = os.path.join(self.temp_dir, 'pkl')
        utils.write_pickle('test', path)
        with gzip.open(path, 'rb') as f:
            self.assertEqual('test', pickle.load(f))

    def test_read_pickle(self):
        path = os.path.join(self.temp_dir, 'pkl')
        with gzip.open(path, 'wb') as f:
            pickle.dump('test', f)
        self.assertEqual('test', utils.read_pickle(path))

    def test_flatten_dict_values(self):
        d = {
            'key1': 'value1',
            'key2': {
                'key3': {
                    'key4': 'value2',
                    'key5': 'value3'
                },
                'key6': 'value4'
            }
        }
        self.assertEqual(['value1', 'value2', 'value3', 'value4'],
                         utils.flatten_dict_values(d))

    def test_flatten_dictionary(self):
        d = {'a': 'b', 'c': 'd', 'e': {'f': 'g', 'h': {'i': 'j'}}}
        self.assertEqual([(('a',), 'b'), (('c',), 'd'), (('e', 'f'), 'g'),
                          (('e', 'h', 'i'), 'j')],
                         list(utils.flatten_dictionary(d)))

    def test_flatten_iter(self):
        lst = [1, [2, 3], [4, [5, 'str']]]
        self.assertEqual([1, 2, 3, 4, 5, 'str'], list(utils.flatten_iter(lst)))

    def test_merge_dictionaries(self):
        d1 = {'a': 'b', 'c': {'d': 'e'}, 'f': 'g'}
        d2 = {'a': 'h', 'c': {'i': 'j'}}
        self.assertEqual({
            'a': 'bh',
            'c': {
                'd': 'eX',
                'i': 'Xj'
            },
            'f': 'gX'
        }, utils.merge_dictionaries(d1, d2, default='X'))
