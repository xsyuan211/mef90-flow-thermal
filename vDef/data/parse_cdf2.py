"""Minimal NetCDF classic / CDF2 reader for Exodus .gen files.

Only the pieces needed for post-processing are implemented: header
(dimensions, global attributes, variable table) + raw data offsets.
All multi-byte values in CDF2 are big-endian.
"""
import struct

NC_DIMENSION = 10
NC_VARIABLE = 11
NC_ATTRIBUTE = 12
NC_BYTE, NC_CHAR, NC_SHORT, NC_INT, NC_FLOAT, NC_DOUBLE = 1, 2, 3, 4, 5, 6
TYPE_SIZE = {1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8}
TYPE_FMT = {1: 'b', 2: 'c', 3: 'h', 4: 'i', 5: 'f', 6: 'd'}


def pad4(n):
    return (n + 3) & ~3


class CDF2:
    def __init__(self, path):
        self.b = open(path, 'rb').read()
        if self.b[:4] != b'CDF\x02':
            raise ValueError('not a CDF2 file: %r' % self.b[:4])
        self.numrecs = struct.unpack_from('>I', self.b, 4)[0]
        self.pos = 8
        self.dims, self.sizes = [], []
        self.vars = []
        self._read_dims()
        self._read_gatts()
        self._read_vars()

    def u32(self):
        v = struct.unpack_from('>I', self.b, self.pos)[0]
        self.pos += 4
        return v

    def u64(self):
        v = struct.unpack_from('>Q', self.b, self.pos)[0]
        self.pos += 8
        return v

    def name(self):
        n = self.u32()
        raw = self.b[self.pos:self.pos + n]
        self.pos += pad4(n)
        return raw.rstrip(b'\x00').decode('ascii', 'replace')

    def _skip_att_list(self):
        tag = self.u32()
        if tag == 0:                      # ABSENT: tag + nelems
            self.u32()
            return []
        if tag != NC_ATTRIBUTE:
            raise ValueError('bad attribute tag %d' % tag)
        n = self.u32()
        out = []
        for _ in range(n):
            nm = self.name()
            typ = self.u32()
            nelems = self.u32()
            sz = pad4(nelems * TYPE_SIZE[typ])
            raw = self.b[self.pos:self.pos + sz]
            self.pos += sz
            out.append((nm, typ, nelems, raw))
        return out

    def _read_dims(self):
        tag = self.u32()
        n = self.u32()
        if tag == 0:
            return
        if tag != NC_DIMENSION:
            raise ValueError('bad dim tag %d' % tag)
        for _ in range(n):
            nm = self.name()
            sz = self.u32()              # 0 = unlimited (record) dimension
            self.dims.append(nm)
            self.sizes.append(sz)

    def _read_gatts(self):
        self.gatts = self._skip_att_list()

    def _read_vars(self):
        tag = self.u32()
        n = self.u32()
        if tag == 0:
            return
        if tag != NC_VARIABLE:
            raise ValueError('bad var tag %d' % tag)
        for _ in range(n):
            nm = self.name()
            ndim = self.u32()
            dimids = [self.u32() for _ in range(ndim)]
            atts = self._skip_att_list()
            typ = self.u32()
            vsize = self.u32()           # CDF2: 32-bit
            begin = self.u64()           # CDF2: 64-bit
            self.vars.append(dict(name=nm, dimids=dimids, type=typ,
                                  vsize=vsize, begin=begin, atts=atts,
                                  is_record=bool(dimids) and self.sizes[dimids[0]] == 0))

    def var(self, name):
        for v in self.vars:
            if v['name'] == name:
                return v
        raise KeyError(name)

    def read_static(self, name):
        v = self.var(name)
        n = v['vsize'] // TYPE_SIZE[v['type']]
        fmt = '>' + TYPE_FMT[v['type']] * n
        return struct.unpack_from(fmt, self.b, v['begin'])

    def read_record(self, name, step):
        v = self.var(name)
        n = v['vsize'] // TYPE_SIZE[v['type']]
        fmt = '>' + TYPE_FMT[v['type']] * n
        # record vars are interleaved: begin points at step 0, stride = RECSZ
        recs = self.record_stride
        return struct.unpack_from(fmt, self.b, v['begin'] + step * recs)

    @property
    def record_stride(self):
        vs = [v for v in self.vars if v['is_record']]
        return sum(pad4(v['vsize']) for v in vs)

    def summary(self):
        print('numrecs', self.numrecs, 'dims', list(zip(self.dims, self.sizes)))
        print('record_stride', self.record_stride,
              'REC0', min(v['begin'] for v in self.vars if v['is_record']))
        print('%-28s %5s %9s %9s %s' % ('name', 'type', 'vsize', 'begin', 'dims'))
        for v in self.vars:
            print('%-28s %5d %9d %9d %s%s' % (
                v['name'], v['type'], v['vsize'], v['begin'],
                [self.dims[d] for d in v['dimids']],
                '  [REC]' if v['is_record'] else ''))


if __name__ == '__main__':
    import sys
    c = CDF2(sys.argv[1] if len(sys.argv) > 1 else
             'D:/SoftBuild/Mef90+Flow+Dy/vDef/data/Beam1x10-tri3_out.gen')
    c.summary()
