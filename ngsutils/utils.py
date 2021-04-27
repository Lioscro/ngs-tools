def is_gzip(path):
    magic = b'\x1f\x8b'
    with open(path, 'rb') as f:
        return magic == f.read(len(magic))
